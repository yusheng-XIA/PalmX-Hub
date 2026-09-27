#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

RUN = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/nrly_hap2_completion_plot_20260812/RUN-OP41-NRLY-HAP2COMP-PLOT-20260812-001")
V2 = Path("${ANALYSIS_DIR}/22_answer_reviews/00_ms/03_V3/01_figure1/00_final_input_resources/05_syri_plotsr_chain12/nrly_contig_rescaffold_20260811/RUN-OP41-NRLY-REFSCAF-20260811-001")
CAND = RUN / "results/03_selected_final/final"
TYPES = {"SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP"}
COLORS = {
    "African_accession": "#506784", "American_oil_palm": "#0072B2",
    "African_oil_palm": "#D55E00", "Phased_dura": "#009E73",
    "Phased_pisifera": "#CC79A7", "Phased_nrly": "#E69F00", "Phased_BK": "#6F58A8",
}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def merge_chr_results(root: Path, output: Path) -> None:
    with output.open("w") as out:
        for chrom in range(1, 17):
            src = root / f"chr{chrom:02d}" / "final/syri.out"
            if not src.is_file() or not src.stat().st_size:
                raise RuntimeError(f"Missing completed SyRI result: {src}")
            with src.open() as inp:
                for line in inp:
                    if line.strip():
                        out.write(line)


def merge_nrly_selected(root: Path, old_source: Path, output: Path) -> None:
    """Use v3 validation except chr04, which was rejected and reverted to v2."""
    with output.open("w") as out:
        for chrom in range(1, 17):
            chrom_id = f"chr{chrom:02d}"
            if chrom_id == "chr04":
                with old_source.open() as inp:
                    for line in inp:
                        if line.startswith("chr04\t"):
                            out.write(line)
                continue
            src = root / chrom_id / "final/syri.out"
            if not src.is_file() or not src.stat().st_size:
                raise RuntimeError(f"Missing completed SyRI result: {src}")
            with src.open() as inp:
                for line in inp:
                    if line.strip():
                        out.write(line)


def filter_syri(src: Path, dst: Path) -> int:
    n = 0
    with src.open() as inp, dst.open("w") as out:
        for line in inp:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 11 and f[10] in TYPES:
                out.write(line); n += 1
    return n


def main() -> None:
    genomes = read_tsv(V2 / "config/Genome_Order.final.tsv")
    genomes = [r for r in genomes if r["Genome_ID"] not in {"EO12", "EG11"}]
    for order, row in enumerate(genomes, 1):
        row["Order"] = str(order)
        if row["Genome_ID"] == "nrly_hap2":
            row["FASTA"] = str(CAND / "nrly_hap2.v3.fa")
    if len(genomes) != 39:
        raise RuntimeError(f"Expected 39 displayed genomes, observed {len(genomes)}")
    write_tsv(RUN / "config/Genome_Order.plot39.tsv", genomes, list(genomes[0]))

    old_pairs = read_tsv(V2 / "config/Pair_Manifest.final.tsv")
    old_by_ids = {(r["Reference_ID"], r["Query_ID"]): r for r in old_pairs}
    final_edges = RUN / "results/06_plot39_final_edges_linked"
    filtered = RUN / "results/07_plot39_plotsr_inputs_linked/filtered_edges"
    final_edges.mkdir(parents=True, exist_ok=True)
    filtered.mkdir(parents=True, exist_ok=True)

    special = {
        ("MZ4_hap2", "American_hap1"): RUN / "results/05_plot39_bridge_edges/edge14_MZ4_hap2__American_hap1",
        ("Africa_hap2", "dura_hap1"): RUN / "results/05_plot39_bridge_edges/edge16_Africa_hap2__dura_hap1",
        ("nrly_hap1", "nrly_hap2"): RUN / "results/04_validation_edges_strict/edge23_nrly_hap1__nrly_hap2",
        ("nrly_hap2", "BK_hap1"): RUN / "results/04_validation_edges_strict/edge24_nrly_hap2__BK_hap1",
    }
    pairs = []
    for edge, (ref, qry) in enumerate(zip(genomes, genomes[1:]), 1):
        ids = (ref["Genome_ID"], qry["Genome_ID"])
        merged = final_edges / f"edge{edge:02d}_{ids[0]}__{ids[1]}.syri.out"
        if ids in special:
            if ids[0] == "nrly_hap1" or ids[0] == "nrly_hap2":
                merge_nrly_selected(special[ids], Path(old_by_ids[ids]["Source_SyRI"]), merged)
            else:
                merge_chr_results(special[ids], merged)
            status = "new_plot39_bridge_20260812" if ids[0] in {"MZ4_hap2", "Africa_hap2"} else "nrly_hap2_v3_strict_20260812"
            dst = filtered / merged.name
            retained = filter_syri(merged, dst)
            if retained == 0:
                raise RuntimeError(f"No plotsr records retained for {ids}")
        elif ids in old_by_ids:
            src = Path(old_by_ids[ids]["Source_SyRI"])
            merged.symlink_to(src)
            old_edge = int(old_by_ids[ids]["Edge"])
            old_filtered = V2 / "results/06_plotsr_inputs/filtered_edges" / f"edge{old_edge:02d}_{ids[0]}__{ids[1]}.syri.out"
            if not old_filtered.is_file():
                raise RuntimeError(f"Missing reviewed filtered edge: {old_filtered}")
            (filtered / merged.name).symlink_to(old_filtered)
            status = old_by_ids[ids]["Status"]
        else:
            raise RuntimeError(f"No SyRI source for adjacent pair {ids}")
        pairs.append({"Edge": edge, "Reference_ID": ids[0], "Query_ID": ids[1], "Status": status, "Source_SyRI": str(merged)})
    write_tsv(RUN / "config/Pair_Manifest.plot39.tsv", pairs, ["Edge", "Reference_ID", "Query_ID", "Status", "Source_SyRI"])

    plot = RUN / "results/07_plot39_plotsr_inputs_linked"
    with (plot / "genomes.txt").open("w") as handle:
        for row in genomes:
            fasta = Path(row["FASTA"])
            fai = Path(str(fasta) + ".fai")
            if not fasta.is_file() or not fai.is_file() or len(fai.read_text().splitlines()) != 16:
                raise RuntimeError(f"Invalid 16-chromosome FASTA/index: {fasta}")
            handle.write(f"{fasta}\t{row['Display_Label']}\tlw:0.85;lc:{COLORS[row['Class']]}\n")
    (plot / "chrord.txt").write_text("".join(f"chr{i:02d}\n" for i in range(1, 17)))
    (plot / "plotsr_pretty.cfg").write_text(
        "legend:T\ngenlegcol:12\nsyncol:#D8DDE3\ninvcol:#E69F00\ntracol:#56B4E9\n"
        "dupcol:#D55E00\nalpha:0.80\nchrmar:0.018\nexmar:0.034\nmarginchr:0.005\n"
        "bbox:0,1.012,0.48,0.30\nbboxmar:0.48\n"
    )
    (plot / "PREPARED").touch()
    print("PASS: 39 genomes, 38 adjacent edges, Mb coordinate mode prepared")


if __name__ == "__main__":
    main()
