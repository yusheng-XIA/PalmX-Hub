#!/usr/bin/env python3
import argparse
import bisect
import csv
import gzip
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path("${ANALYSIS_DIR}/21_MS/06_result/dSVs")
DEFAULT_MINIMAP_DIR = Path("${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/minimap2_results")
DEFAULT_SNP_REF_FAI = Path("${ANALYSIS_DIR}/14_pan_genome/04_SNP_calling/00_renamed_genomes/Africa_hap2.fasta.fai")
DEFAULT_CONSERVED_BED = BASE / "results/archive_old_versions/05_dsv_v3_multi_outgroup/multi_outgroup_conserved_regions/multi_outgroup_conserved_regions.bed"
VALID_BASES = {"A", "C", "G", "T"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
EFFECT_RANK = {
    "Stop_Gained": 100,
    "Stop_Lost": 90,
    "Start_Lost": 80,
    "Missense": 70,
    "Splice_Region": 60,
    "Synonymous": 50,
    "Unknown_CDS_Effect": 40,
    "UTR": 30,
    "Intron": 20,
    "Upstream_2kb": 10,
    "Downstream_2kb": 9,
    "Intergenic": 0,
}
FUNCTIONAL_EFFECTS = {"Missense", "Stop_Gained", "Stop_Lost", "Start_Lost", "Splice_Region", "Conserved_Region_Proxy"}


def parse_args():
    p = argparse.ArgumentParser(description="Reproduce the legacy-compatible dSNP population catalogue for hap38.")
    root = BASE / "results-8.9"
    p.add_argument("--minimap-dir", default=str(DEFAULT_MINIMAP_DIR), help="Fallback directory scan when --sample-manifest is empty.")
    p.add_argument("--sample-manifest", default=str(root / "config/dSNP_Hap38_Input_Manifest.tsv"))
    p.add_argument("--ref-fasta", default=str(BASE / "input/Africa_hap2.fa"))
    p.add_argument("--ref-fai", default=str(BASE / "input/Africa_hap2.fa.fai"))
    p.add_argument("--snp-ref-fai", default=str(DEFAULT_SNP_REF_FAI))
    p.add_argument("--gff3", default=str(BASE / "input/Africa_hap2.EVM.gff3"))
    p.add_argument("--dsv-v5-catalog", default=str(root / "01_core_dsv/sv_catalog.dsv_hap38.tsv"))
    p.add_argument("--conserved-bed", default=str(DEFAULT_CONSERVED_BED))
    p.add_argument("--repeat-bed", default="")
    p.add_argument("--outdir", default=str(root / "05_dSNP_minimap_hap38"))
    p.add_argument("--expected-sample-count", type=int, default=38)
    p.add_argument("--min-qual", type=float, default=30.0)
    p.add_argument("--updown-bp", type=int, default=2000)
    p.add_argument("--splice-bp", type=int, default=2)
    p.add_argument("--sv-breakpoint-window-bp", type=int, default=50)
    p.add_argument("--bin-size", type=int, default=100000)
    p.add_argument("--sort-memory", default="24G")
    p.add_argument("--max-ref-mismatch-records-per-sample", type=int, default=1000,
                   help="Cap mismatch examples only; exact mismatch counts remain uncapped.")
    p.add_argument("--threads", type=int, default=max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))
    p.add_argument("--validate-inputs-only", action="store_true")
    p.add_argument("--allow-existing-outdir", action="store_true")
    return p.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def die(message, code=2):
    print(f"[FATAL] {message}", file=sys.stderr)
    raise SystemExit(code)


def as_int(value, default=0):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value, default=0.0):
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return default


def parse_attrs(text):
    out = {}
    for part in text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def write_rows(path, fieldnames, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "NA") for key in fieldnames})


def discover_samples(minimap_dir):
    rows = []
    for sample_dir in sorted(Path(minimap_dir).iterdir()):
        if not sample_dir.is_dir() or sample_dir.name in {"scripts", "logs"}:
            continue
        sample = sample_dir.name
        snp = sample_dir / f"{sample}_vs_Africa_hap2.snp.txt"
        summary = sample_dir / f"{sample}_vs_Africa_hap2.summary.txt"
        rows.append({"Sample": sample, "SNP_File": str(snp), "Summary_File": str(summary)})
    return rows


def load_sample_manifest(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        die(f"Sample manifest missing or empty: {path}", code=3)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sample", "SNP_File", "Summary_File"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            die(f"Sample manifest missing columns: {sorted(missing)}", code=3)
        rows = [{"Sample": row["Sample"], "SNP_File": row["SNP_File"], "Summary_File": row["Summary_File"]} for row in reader]
    names = [row["Sample"] for row in rows]
    if len(names) != len(set(names)):
        die("Duplicate sample names in sample manifest", code=3)
    return rows


def parse_summary(path):
    metrics = {}
    with open(path) as handle:
        for line in handle:
            text = line.strip()
            if not text or ":" not in text:
                continue
            key, val = text.split(":", 1)
            metrics[key.strip()] = val.strip()
    return {
        "Summary_Total": as_int(metrics.get("Total"), "NA"),
        "Summary_SNP": as_int(metrics.get("SNP"), "NA"),
        "Summary_InDel": as_int(metrics.get("InDel"), "NA"),
    }


def read_fai(path):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def build_chrom_map(snp_lengths, ref_lengths):
    ref_by_length = defaultdict(list)
    for chrom, length in ref_lengths.items():
        ref_by_length[length].append(chrom)
    chrom_map = {}
    ambiguous = []
    missing = []
    for chrom, length in snp_lengths.items():
        if chrom in ref_lengths and ref_lengths[chrom] == length:
            chrom_map[chrom] = chrom
            continue
        matches = sorted(ref_by_length.get(length, []))
        if len(matches) == 1:
            chrom_map[chrom] = matches[0]
        elif len(matches) > 1:
            ambiguous.append((chrom, length, ",".join(matches)))
        else:
            missing.append((chrom, length))
    if ambiguous or missing:
        details = []
        if ambiguous:
            details.append(f"ambiguous={ambiguous[:5]}")
        if missing:
            details.append(f"missing={missing[:5]}")
        die("Cannot build unambiguous SNP-to-reference chromosome map: " + " | ".join(details), code=4)
    return chrom_map


def write_chrom_map(path, chrom_map, snp_lengths, ref_lengths):
    rows = []
    for source in sorted(chrom_map):
        target = chrom_map[source]
        rows.append({
            "SNP_Chrom": source,
            "Reference_Chrom": target,
            "SNP_Length": snp_lengths[source],
            "Reference_Length": ref_lengths[target],
            "Status": "PASS" if snp_lengths[source] == ref_lengths[target] else "FAIL",
        })
    write_rows(path, ["SNP_Chrom", "Reference_Chrom", "SNP_Length", "Reference_Length", "Status"], rows)


def read_fasta(path, chrom_whitelist):
    seqs = {}
    chrom = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if chrom is not None and chrom in chrom_whitelist:
                    seqs[chrom] = "".join(chunks).upper()
                chrom = line[1:].split()[0]
                chunks = []
            elif chrom in chrom_whitelist:
                chunks.append(line.strip())
        if chrom is not None and chrom in chrom_whitelist:
            seqs[chrom] = "".join(chunks).upper()
    return seqs


def write_n_regions(seqs, out_path):
    rows = 0
    with open(out_path, "w") as out:
        for chrom, seq in seqs.items():
            start = None
            for i, base in enumerate(seq):
                if base == "N" and start is None:
                    start = i
                elif base != "N" and start is not None:
                    out.write(f"{chrom}\t{start}\t{i}\n")
                    rows += 1
                    start = None
            if start is not None:
                out.write(f"{chrom}\t{start}\t{len(seq)}\n")
                rows += 1
    return rows


class IntervalIndex:
    def __init__(self, intervals=None, bin_size=100000):
        self.bin_size = bin_size
        self.intervals = defaultdict(list)
        if intervals:
            for chrom, start, end, payload in intervals:
                if end > start:
                    self.intervals[chrom].append((int(start), int(end), payload))
        self.bins = defaultdict(lambda: defaultdict(list))
        for chrom, chrom_intervals in self.intervals.items():
            chrom_intervals.sort(key=lambda x: (x[0], x[1]))
            for idx, (start, end, payload) in enumerate(chrom_intervals):
                for b in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
                    self.bins[chrom][b].append(idx)

    def hits(self, chrom, pos0):
        out = []
        b = pos0 // self.bin_size
        for idx in self.bins.get(chrom, {}).get(b, []):
            start, end, payload = self.intervals[chrom][idx]
            if start <= pos0 < end:
                out.append((start, end, payload))
        return out


def read_bed_index(path, bin_size, label="BED"):
    if not path:
        return IntervalIndex(bin_size=bin_size)
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return IntervalIndex(bin_size=bin_size)
    intervals = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            intervals.append((fields[0], as_int(fields[1]), as_int(fields[2]), fields[3] if len(fields) > 3 else label))
    return IntervalIndex(intervals, bin_size=bin_size)


def read_sv_breakpoint_index(path, window, bin_size):
    intervals = []
    if not Path(path).exists():
        return IntervalIndex(bin_size=bin_size)
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row.get("Chrom") or row.get("Original_Chrom")
            if not chrom or chrom == "NA":
                continue
            start0 = as_int(row.get("Interval_Start0", row.get("Start", row.get("Pos", 0))), None)
            end0 = as_int(row.get("Interval_End0", row.get("End", row.get("Pos", 0))), None)
            if start0 is None or end0 is None:
                continue
            key = row.get("SV_Key", row.get("SV_ID", "SV"))
            intervals.append((chrom, max(0, start0 - window), start0 + window + 1, key))
            intervals.append((chrom, max(0, end0 - window), end0 + window + 1, key))
    return IntervalIndex(intervals, bin_size=bin_size)


class Annotation:
    def __init__(self, gff3, seqs, updown_bp, splice_bp, bin_size):
        self.seqs = seqs
        self.updown_bp = updown_bp
        self.splice_bp = splice_bp
        self.genes = []
        self.exons = []
        self.cds = []
        self.upstream = []
        self.downstream = []
        self.mrna_to_gene = {}
        self.tx_strand = {}
        self.tx_chrom = {}
        self.tx_cds = defaultdict(list)
        self.tx_exons = defaultdict(list)
        self.cds_cache = {}
        self._read_gff(gff3)
        self.gene_index = IntervalIndex(self.genes, bin_size)
        self.exon_index = IntervalIndex(self.exons, bin_size)
        self.cds_index = IntervalIndex(self.cds, bin_size)
        self.upstream_index = IntervalIndex(self.upstream, bin_size)
        self.downstream_index = IntervalIndex(self.downstream, bin_size)

    def _read_gff(self, path):
        gene_bounds = {}
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _source, feature, start, end, _score, strand, phase, attrs_text = fields
                start0 = int(start) - 1
                end0 = int(end)
                attrs = parse_attrs(attrs_text)
                if feature == "gene":
                    gene_id = attrs.get("ID", attrs.get("Name", "NA"))
                    gene_bounds[gene_id] = (chrom, start0, end0, strand)
                    self.genes.append((chrom, start0, end0, {"Gene_ID": gene_id, "Strand": strand}))
                elif feature in {"mRNA", "transcript"}:
                    tx = attrs.get("ID", "NA")
                    gene_id = attrs.get("Parent", "NA").split(",")[0]
                    self.mrna_to_gene[tx] = gene_id
                    self.tx_strand[tx] = strand
                    self.tx_chrom[tx] = chrom
                elif feature == "exon":
                    for tx in attrs.get("Parent", "NA").split(","):
                        gene_id = self.mrna_to_gene.get(tx, "NA")
                        payload = {"Transcript_ID": tx, "Gene_ID": gene_id, "Strand": strand}
                        self.exons.append((chrom, start0, end0, payload))
                        self.tx_exons[tx].append((start0, end0))
                elif feature == "CDS":
                    for tx in attrs.get("Parent", "NA").split(","):
                        gene_id = self.mrna_to_gene.get(tx, "NA")
                        payload = {"Transcript_ID": tx, "Gene_ID": gene_id, "Strand": strand, "Phase": phase}
                        self.cds.append((chrom, start0, end0, payload))
                        self.tx_cds[tx].append((chrom, start0, end0, strand))
                        self.tx_chrom[tx] = chrom
                        self.tx_strand[tx] = strand
        for gene_id, (chrom, start0, end0, strand) in gene_bounds.items():
            if strand == "-":
                self.upstream.append((chrom, end0, end0 + self.updown_bp, {"Gene_ID": gene_id, "Strand": strand}))
                self.downstream.append((chrom, max(0, start0 - self.updown_bp), start0, {"Gene_ID": gene_id, "Strand": strand}))
            else:
                self.upstream.append((chrom, max(0, start0 - self.updown_bp), start0, {"Gene_ID": gene_id, "Strand": strand}))
                self.downstream.append((chrom, end0, end0 + self.updown_bp, {"Gene_ID": gene_id, "Strand": strand}))

    def _tx_model(self, tx):
        if tx in self.cds_cache:
            return self.cds_cache[tx]
        parts = self.tx_cds.get(tx, [])
        if not parts:
            self.cds_cache[tx] = None
            return None
        chrom = parts[0][0]
        strand = parts[0][3]
        ordered = sorted(parts, key=lambda x: x[1], reverse=(strand == "-"))
        seq_parts = []
        genomic_positions = []
        for _chrom, start, end, _strand in ordered:
            piece = self.seqs.get(chrom, "")[start:end]
            positions = list(range(start, end))
            if strand == "-":
                piece = piece.translate(RC)[::-1].upper()
                positions = positions[::-1]
            seq_parts.append(piece)
            genomic_positions.extend(positions)
        seq = "".join(seq_parts).upper()
        pos_to_cds = {pos: i for i, pos in enumerate(genomic_positions)}
        model = {"Chrom": chrom, "Strand": strand, "Seq": seq, "Pos_To_CDS": pos_to_cds}
        self.cds_cache[tx] = model
        return model

    def _cds_effect(self, chrom, pos0, ref, alt):
        calls = []
        for _start, _end, payload in self.cds_index.hits(chrom, pos0):
            tx = payload["Transcript_ID"]
            model = self._tx_model(tx)
            if not model:
                continue
            cds_i = model["Pos_To_CDS"].get(pos0)
            if cds_i is None:
                continue
            codon_start = (cds_i // 3) * 3
            codon = model["Seq"][codon_start:codon_start + 3]
            if len(codon) != 3 or any(base not in VALID_BASES for base in codon):
                calls.append(("Unknown_CDS_Effect", payload["Gene_ID"], tx, "NA", "NA", "NA"))
                continue
            alt_base = alt if model["Strand"] != "-" else alt.translate(RC)
            ref_base = ref if model["Strand"] != "-" else ref.translate(RC)
            rel = cds_i - codon_start
            if codon[rel] != ref_base:
                calls.append(("Unknown_CDS_Effect", payload["Gene_ID"], tx, codon, "NA", "NA"))
                continue
            alt_codon = codon[:rel] + alt_base + codon[rel + 1:]
            ref_aa = CODON_TABLE.get(codon, "X")
            alt_aa = CODON_TABLE.get(alt_codon, "X")
            if ref_aa != "*" and alt_aa == "*":
                effect = "Stop_Gained"
            elif ref_aa == "*" and alt_aa != "*":
                effect = "Stop_Lost"
            elif codon_start == 0 and ref_aa == "M" and alt_aa != "M":
                effect = "Start_Lost"
            elif ref_aa == alt_aa:
                effect = "Synonymous"
            elif ref_aa != alt_aa:
                effect = "Missense"
            else:
                effect = "Unknown_CDS_Effect"
            aa_change = f"{ref_aa}{cds_i // 3 + 1}{alt_aa}"
            calls.append((effect, payload["Gene_ID"], tx, codon, alt_codon, aa_change))
        if not calls:
            return None
        calls.sort(key=lambda x: EFFECT_RANK.get(x[0], -1), reverse=True)
        return calls[0]

    def _splice_hit(self, chrom, pos0):
        for start, end, payload in self.exon_index.hits(chrom, pos0):
            if min(abs(pos0 - start), abs(pos0 - (end - 1))) <= self.splice_bp:
                return payload
        return None

    def annotate(self, chrom, pos, ref, alt, conserved_proxy):
        pos0 = int(pos) - 1
        cds = self._cds_effect(chrom, pos0, ref, alt)
        splice = self._splice_hit(chrom, pos0)
        gene_hits = self.gene_index.hits(chrom, pos0)
        exon_hits = self.exon_index.hits(chrom, pos0)
        upstream_hits = self.upstream_index.hits(chrom, pos0)
        downstream_hits = self.downstream_index.hits(chrom, pos0)
        gene_ids = sorted({h[2].get("Gene_ID", "NA") for h in gene_hits + upstream_hits + downstream_hits + exon_hits})
        tx_ids = sorted({h[2].get("Transcript_ID", "NA") for h in exon_hits if h[2].get("Transcript_ID", "NA") != "NA"})
        if cds:
            effect, gene_id, tx, ref_codon, alt_codon, aa_change = cds
            context = "CDS"
            gene_ids = [gene_id]
            tx_ids = [tx]
        elif splice:
            effect = "Splice_Region"
            context = "Splice_Region"
            ref_codon = alt_codon = aa_change = "NA"
        elif exon_hits:
            effect = "UTR"
            context = "UTR"
            ref_codon = alt_codon = aa_change = "NA"
        elif gene_hits:
            effect = "Intron"
            context = "Intron"
            ref_codon = alt_codon = aa_change = "NA"
        elif upstream_hits:
            effect = "Upstream_2kb"
            context = "Upstream_2kb"
            ref_codon = alt_codon = aa_change = "NA"
        elif downstream_hits:
            effect = "Downstream_2kb"
            context = "Downstream_2kb"
            ref_codon = alt_codon = aa_change = "NA"
        else:
            effect = "Intergenic"
            context = "Intergenic"
            ref_codon = alt_codon = aa_change = "NA"
        if conserved_proxy and effect not in {"Stop_Gained", "Stop_Lost", "Start_Lost", "Missense", "Splice_Region"}:
            functional_class = "Conserved_Region_Proxy"
        else:
            functional_class = effect
        return {
            "Gene_Context": context,
            "Functional_Effect": effect,
            "Functional_Class": functional_class,
            "Gene_IDs": ";".join(gene_ids) if gene_ids else "NA",
            "Transcript_IDs": ";".join(tx_ids) if tx_ids else "NA",
            "Ref_Codon": ref_codon,
            "Alt_Codon": alt_codon,
            "AA_Change": aa_change,
            "Splice_Region": "Yes" if effect == "Splice_Region" else "No",
        }


def validate_inputs(args, samples, outdir):
    rows = []
    failures = []
    if len(samples) != args.expected_sample_count:
        failures.append(f"Expected {args.expected_sample_count} sample dirs, found {len(samples)}")
    for row in samples:
        sample = row["Sample"]
        snp = Path(row["SNP_File"])
        summary = Path(row["Summary_File"])
        status = "PASS"
        notes = []
        if not snp.exists() or snp.stat().st_size == 0:
            status = "FAIL"
            notes.append("missing_or_empty_snp")
        if not summary.exists() or summary.stat().st_size == 0:
            status = "FAIL"
            notes.append("missing_or_empty_summary")
        metrics = parse_summary(summary) if summary.exists() and summary.stat().st_size > 0 else {}
        rows.append({
            **row,
            "SNP_File_Size": snp.stat().st_size if snp.exists() else 0,
            "Summary_File_Size": summary.stat().st_size if summary.exists() else 0,
            "Summary_Total": metrics.get("Summary_Total", "NA"),
            "Summary_SNP": metrics.get("Summary_SNP", "NA"),
            "Summary_InDel": metrics.get("Summary_InDel", "NA"),
            "Status": status,
            "Notes": ";".join(notes) if notes else "OK",
        })
        if status != "PASS":
            failures.append(f"{sample}: {rows[-1]['Notes']}")
    for path in [args.ref_fasta, args.ref_fai, args.snp_ref_fai, args.gff3, args.dsv_v5_catalog]:
        if not Path(path).exists() or Path(path).stat().st_size == 0:
            failures.append(f"Required input missing or empty: {path}")
    if args.conserved_bed and (not Path(args.conserved_bed).exists() or Path(args.conserved_bed).stat().st_size == 0):
        print(f"[WARN] conserved BED not found or empty; Conserved_Region_Proxy will be No: {args.conserved_bed}", file=sys.stderr)
    if outdir:
        write_rows(outdir / "sample_manifest.tsv", [
            "Sample", "SNP_File", "Summary_File", "SNP_File_Size", "Summary_File_Size",
            "Summary_Total", "Summary_SNP", "Summary_InDel", "Status", "Notes",
        ], rows)
    if failures:
        die("Input validation failed: " + " | ".join(failures), code=3)
    return rows


def normalize_records(args, samples, seqs, lengths, chrom_map, tmp_records, qa_dir):
    audit_rows = []
    ref_conflicts_pre = qa_dir / "ref_mismatch_records.tsv"
    with open(tmp_records, "w") as out, open(ref_conflicts_pre, "w") as ref_bad:
        out.write("Chrom\tPos\tRef\tAlt\tSample\tQual\tFilter\n")
        ref_bad.write("Sample\tSNP_Chrom\tReference_Chrom\tPos\tInput_Ref\tReference_Base\tAlt\tQual\tFilter\n")
        for sample_row in samples:
            sample = sample_row["Sample"]
            snp_path = sample_row["SNP_File"]
            seen = set()
            counts = Counter()
            with open_text(snp_path) as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    counts["Raw_SNP_Rows"] += 1
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 8:
                        counts["Malformed_Rows"] += 1
                        continue
                    source_chrom, pos_text, _id, ref, alt, qual_text, filt, _info = fields[:8]
                    chrom = chrom_map.get(source_chrom)
                    if not chrom:
                        counts["Unmapped_Chrom"] += 1
                        continue
                    if chrom != source_chrom:
                        counts["Chrom_Mapped"] += 1
                    ref = ref.upper()
                    alt = alt.upper()
                    pos = as_int(pos_text, None)
                    qual = as_float(qual_text, 0.0)
                    if pos is None or chrom not in lengths or pos < 1 or pos > lengths[chrom]:
                        counts["Out_Of_Bounds"] += 1
                        continue
                    if len(ref) != 1 or len(alt) != 1 or ref not in VALID_BASES or alt not in VALID_BASES:
                        counts["Non_SNV_or_Non_ACGT"] += 1
                        continue
                    ref_base = seqs[chrom][pos - 1]
                    if ref_base == "N":
                        counts["Reference_N"] += 1
                        continue
                    if ref_base != ref:
                        counts["Reference_Mismatch"] += 1
                        if counts["Reference_Mismatch_Records_Written"] < args.max_ref_mismatch_records_per_sample:
                            ref_bad.write(f"{sample}\t{source_chrom}\t{chrom}\t{pos}\t{ref}\t{ref_base}\t{alt}\t{qual_text}\t{filt}\n")
                            counts["Reference_Mismatch_Records_Written"] += 1
                        continue
                    key = f"{chrom}\t{pos}\t{ref}\t{alt}"
                    if key in seen:
                        counts["Duplicate_Within_Sample"] += 1
                        continue
                    seen.add(key)
                    counts["Valid_SNV_RefChecked"] += 1
                    out.write(f"{key}\t{sample}\t{qual:.6f}\t{filt}\n")
            metrics = parse_summary(sample_row["Summary_File"])
            audit_rows.append({
                "Sample": sample,
                "Summary_SNP": metrics.get("Summary_SNP", "NA"),
                "Raw_SNP_Rows": counts["Raw_SNP_Rows"],
                "Valid_SNV_RefChecked": counts["Valid_SNV_RefChecked"],
                "Non_SNV_or_Non_ACGT": counts["Non_SNV_or_Non_ACGT"],
                "Out_Of_Bounds": counts["Out_Of_Bounds"],
                "Reference_N": counts["Reference_N"],
                "Reference_Mismatch": counts["Reference_Mismatch"],
                "Reference_Mismatch_Records_Written": counts["Reference_Mismatch_Records_Written"],
                "Chrom_Mapped": counts["Chrom_Mapped"],
                "Unmapped_Chrom": counts["Unmapped_Chrom"],
                "Duplicate_Within_Sample": counts["Duplicate_Within_Sample"],
                "Malformed_Rows": counts["Malformed_Rows"],
                "Raw_vs_Summary_Delta": counts["Raw_SNP_Rows"] - as_int(metrics.get("Summary_SNP"), counts["Raw_SNP_Rows"]),
            })
            print(f"[INFO] normalized {sample}: raw={counts['Raw_SNP_Rows']} valid={counts['Valid_SNV_RefChecked']}", flush=True)
    write_rows(qa_dir / "input_count_audit.tsv", [
        "Sample", "Summary_SNP", "Raw_SNP_Rows", "Valid_SNV_RefChecked",
        "Non_SNV_or_Non_ACGT", "Out_Of_Bounds", "Reference_N", "Reference_Mismatch", "Reference_Mismatch_Records_Written",
        "Chrom_Mapped", "Unmapped_Chrom", "Duplicate_Within_Sample", "Malformed_Rows", "Raw_vs_Summary_Delta",
    ], audit_rows)
    return audit_rows


def sort_records(args, tmp_records, sorted_records, tmpdir):
    cmd = [
        "sort",
        "-T", str(tmpdir),
        "--parallel", str(max(1, args.threads)),
        "-S", args.sort_memory,
        "-k1,1",
        "-k2,2n",
        "-k3,3",
        "-k4,4",
        str(tmp_records),
    ]
    print("[INFO] running sort: " + " ".join(cmd), flush=True)
    with open(sorted_records, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def iter_position_groups(sorted_records):
    with open(sorted_records) as handle:
        header = next(handle, None)
        current_pos = None
        groups = defaultdict(list)
        for line in handle:
            chrom, pos, ref, alt, sample, qual, filt = line.rstrip("\n").split("\t")
            pos_key = (chrom, pos)
            if current_pos is None:
                current_pos = pos_key
            if pos_key != current_pos:
                yield current_pos, groups
                current_pos = pos_key
                groups = defaultdict(list)
            groups[(chrom, pos, ref, alt)].append((sample, as_float(qual), filt))
        if current_pos is not None:
            yield current_pos, groups


def emit_outputs(args, sorted_records, outdir, qa_dir, annotation, sv_index, conserved_index, repeat_index, lengths):
    catalog_fields = [
        "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Sample_Count", "Frequency", "Samples",
        "Mean_Qual", "Filter", "Pass_Count", "Source", "Near_SV_Breakpoint",
        "SV_Breakpoint_Keys", "Repeat_Overlap", "Conserved_Region_Proxy",
    ]
    ann_fields = [
        "SNP_ID", "Chrom", "Pos", "Ref", "Alt", "Sample_Count", "Frequency", "Samples",
        "Mean_Qual", "Filter", "Gene_Context", "Functional_Effect", "Functional_Class",
        "Splice_Region", "Gene_IDs", "Transcript_IDs", "Ref_Codon", "Alt_Codon", "AA_Change",
        "Near_SV_Breakpoint", "Repeat_Overlap", "Conserved_Region_Proxy", "Polarity", "Polarity_Evidence",
    ]
    dsnp_fields = ann_fields + ["Candidate_Class", "dSNP_v1_Flag", "dSNP_v1_Evidence"]
    summary_counts = Counter()
    sample_burden = Counter()
    chrom_counts = Counter()
    effect_counts = Counter()
    conflict_count = 0
    snp_index = 0
    rare_index = 0
    functional_index = 0
    dsnp_index = 0
    with open(outdir / "snp_population_catalog.tsv", "w", newline="") as catalog, \
        open(outdir / "rare_snp_candidates.tsv", "w", newline="") as rare, \
        open(outdir / "snp_functional_annotation.tsv", "w", newline="") as ann, \
        open(outdir / "rare_functional_snp.tsv", "w", newline="") as rare_func, \
        open(outdir / "dsnp_v1_candidates.tsv", "w", newline="") as dsnp, \
        open(qa_dir / "ref_conflict_sites.tsv", "w", newline="") as conflicts:
        catalog_w = csv.DictWriter(catalog, fieldnames=catalog_fields, delimiter="\t", lineterminator="\n")
        rare_w = csv.DictWriter(rare, fieldnames=catalog_fields, delimiter="\t", lineterminator="\n")
        ann_w = csv.DictWriter(ann, fieldnames=ann_fields, delimiter="\t", lineterminator="\n")
        rare_func_w = csv.DictWriter(rare_func, fieldnames=dsnp_fields, delimiter="\t", lineterminator="\n")
        dsnp_w = csv.DictWriter(dsnp, fieldnames=dsnp_fields, delimiter="\t", lineterminator="\n")
        conflict_w = csv.DictWriter(conflicts, fieldnames=["Chrom", "Pos", "Refs", "Alt_Count", "Samples"], delimiter="\t", lineterminator="\n")
        for writer in [catalog_w, rare_w, ann_w, rare_func_w, dsnp_w, conflict_w]:
            writer.writeheader()
        for (chrom, pos), groups in iter_position_groups(sorted_records):
            refs = sorted({key[2] for key in groups})
            if len(refs) > 1:
                conflict_count += 1
                samples = sorted({rec[0] for recs in groups.values() for rec in recs})
                conflict_w.writerow({"Chrom": chrom, "Pos": pos, "Refs": ";".join(refs), "Alt_Count": len(groups), "Samples": ";".join(samples)})
                continue
            for key, recs in sorted(groups.items(), key=lambda item: item[0]):
                chrom, pos, ref, alt = key
                samples = sorted({r[0] for r in recs})
                sample_count = len(samples)
                quals = [r[1] for r in recs]
                filters = sorted({r[2] for r in recs})
                pass_count = sum(1 for r in recs if r[2] == "PASS")
                pos0 = int(pos) - 1
                sv_hits = sv_index.hits(chrom, pos0)
                conserved = bool(conserved_index.hits(chrom, pos0))
                repeat = bool(repeat_index.hits(chrom, pos0))
                snp_index += 1
                snp_id = f"SNP{snp_index:012d}"
                row = {
                    "SNP_ID": snp_id,
                    "Chrom": chrom,
                    "Pos": pos,
                    "Ref": ref,
                    "Alt": alt,
                    "Sample_Count": sample_count,
                    "Frequency": f"{sample_count / args.expected_sample_count:.8f}",
                    "Samples": ";".join(samples),
                    "Mean_Qual": f"{sum(quals) / len(quals):.4f}" if quals else "0.0000",
                    "Filter": "PASS" if len(filters) == 1 and filters[0] == "PASS" else ";".join(filters),
                    "Pass_Count": pass_count,
                    "Source": "minimap2_paftools",
                    "Near_SV_Breakpoint": "Yes" if sv_hits else "No",
                    "SV_Breakpoint_Keys": ";".join(sorted({h[2] for h in sv_hits})) if sv_hits else "NA",
                    "Repeat_Overlap": "Yes" if repeat else "No",
                    "Conserved_Region_Proxy": "Yes" if conserved else "No",
                }
                catalog_w.writerow(row)
                summary_counts["Total_SNV_Catalog"] += 1
                chrom_counts[(chrom, "Catalog")] += 1
                if sample_count > args.expected_sample_count:
                    summary_counts["Sample_Count_Overflow"] += 1
                if sample_count == 1:
                    rare_index += 1
                    rare_w.writerow(row)
                    summary_counts["Rare_Singleton_SNV"] += 1
                    sample_burden[(samples[0], "Rare_Singleton_SNV")] += 1
                    ann_row = dict(row)
                    ann_row.update(annotation.annotate(chrom, int(pos), ref, alt, conserved))
                    ann_row["Polarity"] = "Unknown_Polarity"
                    ann_row["Polarity_Evidence"] = "No_outgroup_base_polarization_in_minimap_v1"
                    ann_w.writerow({field: ann_row.get(field, "NA") for field in ann_fields})
                    effect_counts[ann_row["Functional_Class"]] += 1
                    is_basic_pass = row["Filter"] == "PASS" and as_float(row["Mean_Qual"]) >= args.min_qual and row["Near_SV_Breakpoint"] == "No"
                    functional_evidence = ann_row["Functional_Class"] in FUNCTIONAL_EFFECTS
                    if is_basic_pass and functional_evidence:
                        functional_index += 1
                        rf = dict(ann_row)
                        rf["Candidate_Class"] = "rare_functional_snp"
                        rf["dSNP_v1_Flag"] = "No"
                        rf["dSNP_v1_Evidence"] = "Rare_functional_without_ALT_Derived_or_conserved_proxy"
                        if conserved:
                            rf["dSNP_v1_Flag"] = "Yes"
                            rf["dSNP_v1_Evidence"] = "Conserved_Region_Proxy"
                        rare_func_w.writerow({field: rf.get(field, "NA") for field in dsnp_fields})
                        sample_burden[(samples[0], "Rare_Functional_SNV")] += 1
                        chrom_counts[(chrom, "Rare_Functional_SNV")] += 1
                        if rf["dSNP_v1_Flag"] == "Yes":
                            dsnp_index += 1
                            dsnp_w.writerow({field: rf.get(field, "NA") for field in dsnp_fields})
                            sample_burden[(samples[0], "Putative_dSNP_v1")] += 1
                            chrom_counts[(chrom, "Putative_dSNP_v1")] += 1
    summary_counts["Ref_Conflict_Sites"] = conflict_count
    summary_counts["Rare_Functional_SNV"] = functional_index
    summary_counts["Putative_dSNP_v1"] = dsnp_index
    return summary_counts, sample_burden, chrom_counts, effect_counts


def write_summary(outdir, summary_counts, sample_burden, chrom_counts, effect_counts, lengths, args):
    rows = []
    for metric, value in sorted(summary_counts.items()):
        rows.append({"Summary_Type": "Global", "Group": metric, "Count": value, "Value": value, "Notes": "NA"})
    for (sample, klass), value in sorted(sample_burden.items()):
        rows.append({"Summary_Type": "Sample_Burden", "Group": f"{sample}|{klass}", "Count": value, "Value": value, "Notes": "NA"})
    for (chrom, klass), value in sorted(chrom_counts.items()):
        length = lengths.get(chrom, 0)
        density = f"{value / (length / 1000000.0):.6f}" if length else "NA"
        rows.append({"Summary_Type": "Chromosome_Density", "Group": f"{chrom}|{klass}", "Count": value, "Value": density, "Notes": "Value is count per Mb"})
    for effect, value in sorted(effect_counts.items()):
        rows.append({"Summary_Type": "Functional_Class", "Group": effect, "Count": value, "Value": value, "Notes": "Singleton_annotations_only"})
    rows.append({"Summary_Type": "Method", "Group": "Polarity", "Count": "NA", "Value": "Unknown_Polarity", "Notes": "No outgroup base polarization was applied in v1"})
    rows.append({"Summary_Type": "Method", "Group": "Rare_Definition", "Count": "NA", "Value": "Sample_Count == 1", "Notes": f"Frequency=1/{args.expected_sample_count}"})
    write_rows(outdir / "dsnp_v1_summary.tsv", ["Summary_Type", "Group", "Count", "Value", "Notes"], rows)


def write_integrity(outdir, qa_dir, summary_counts, samples, input_audit, args):
    rows = []
    def add(check, status, observed, expected, notes="NA"):
        rows.append({"Check": check, "Status": status, "Observed": observed, "Expected": expected, "Notes": notes})
    add("sample_count", "PASS" if len(samples) == args.expected_sample_count else "FAIL", len(samples), args.expected_sample_count)
    add("sample_count_overflow", "PASS" if summary_counts.get("Sample_Count_Overflow", 0) == 0 else "FAIL", summary_counts.get("Sample_Count_Overflow", 0), 0)
    add("ref_conflict_sites", "WARN" if summary_counts.get("Ref_Conflict_Sites", 0) else "PASS", summary_counts.get("Ref_Conflict_Sites", 0), 0, "Conflicted positions are excluded from main outputs")
    add("total_catalog_nonempty", "PASS" if summary_counts.get("Total_SNV_Catalog", 0) > 0 else "FAIL", summary_counts.get("Total_SNV_Catalog", 0), ">0")
    add("polarity_caveat", "WARN", "Unknown_Polarity", "ALT_Derived evidence absent", "Strict dSNP set is limited to conserved proxy in v1")
    raw_total = sum(as_int(row.get("Raw_SNP_Rows")) for row in input_audit)
    valid_total = sum(as_int(row.get("Valid_SNV_RefChecked")) for row in input_audit)
    mismatch_total = sum(as_int(row.get("Reference_Mismatch")) for row in input_audit)
    mismatch_written = sum(as_int(row.get("Reference_Mismatch_Records_Written")) for row in input_audit)
    add("legacy_reference_match_fraction", "WARN", f"{valid_total}/{raw_total}", "Legacy-compatible rule",
        f"Reference mismatches excluded={mismatch_total}; see reports/dSNP_Reference_Allele_Audit.md")
    add("reference_mismatch_example_cap", "PASS" if mismatch_written <= len(samples) * args.max_ref_mismatch_records_per_sample else "FAIL",
        mismatch_written, f"<={len(samples) * args.max_ref_mismatch_records_per_sample}", "Counts are exact; only example rows are capped")
    write_rows(qa_dir / "final_integrity_summary.tsv", ["Check", "Status", "Observed", "Expected", "Notes"], rows)
    if any(row["Status"] == "FAIL" for row in rows):
        die("Final integrity checks contain FAIL rows", code=5)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    qa_dir = outdir / "qa"
    tmpdir = outdir / "tmp"
    if outdir.exists() and any(outdir.iterdir()) and not args.allow_existing_outdir and not args.validate_inputs_only:
        die(f"Output directory exists and is non-empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    qa_dir.mkdir(parents=True, exist_ok=True)
    if not args.validate_inputs_only:
        tmpdir.mkdir(parents=True, exist_ok=True)
    samples = load_sample_manifest(args.sample_manifest) if args.sample_manifest else discover_samples(args.minimap_dir)
    manifest_rows = validate_inputs(args, samples, outdir)
    if args.validate_inputs_only:
        print(f"[INFO] input validation OK for {len(samples)} samples")
        return
    lengths = read_fai(args.ref_fai)
    snp_lengths = read_fai(args.snp_ref_fai)
    chrom_map = build_chrom_map(snp_lengths, lengths)
    write_chrom_map(qa_dir / "chrom_name_map.tsv", chrom_map, snp_lengths, lengths)
    print(f"[INFO] chromosome map ready for {len(chrom_map)} SNP reference sequences", flush=True)
    print("[INFO] loading reference FASTA into memory for base/N checks", flush=True)
    seqs = read_fasta(args.ref_fasta, set(lengths))
    missing = sorted(set(lengths) - set(seqs))
    if missing:
        die(f"Reference FASTA missing sequences listed in FAI: {missing[:5]}")
    n_count = write_n_regions(seqs, qa_dir / "reference_n_regions.bed")
    print(f"[INFO] wrote N-region BED with {n_count} intervals", flush=True)
    tmp_records = tmpdir / "all_snp_records.normalized.tsv"
    sorted_records = tmpdir / "all_snp_records.sorted.tsv"
    input_audit = normalize_records(args, manifest_rows, seqs, lengths, chrom_map, tmp_records, qa_dir)
    sort_records(args, tmp_records, sorted_records, tmpdir)
    print("[INFO] loading annotation and interval indexes", flush=True)
    annotation = Annotation(args.gff3, seqs, args.updown_bp, args.splice_bp, args.bin_size)
    sv_index = read_sv_breakpoint_index(args.dsv_v5_catalog, args.sv_breakpoint_window_bp, args.bin_size)
    conserved_index = read_bed_index(args.conserved_bed, args.bin_size, label="Conserved_Region")
    repeat_index = read_bed_index(args.repeat_bed, args.bin_size, label="Repeat")
    summary_counts, sample_burden, chrom_counts, effect_counts = emit_outputs(
        args, sorted_records, outdir, qa_dir, annotation, sv_index, conserved_index, repeat_index, lengths
    )
    write_summary(outdir, summary_counts, sample_burden, chrom_counts, effect_counts, lengths, args)
    write_integrity(outdir, qa_dir, summary_counts, samples, input_audit, args)
    shutil.rmtree(tmpdir)
    (outdir / "_DONE").write_text("dSNP minimap v1 completed\n")
    print(f"[INFO] done: {outdir}")


if __name__ == "__main__":
    main()
