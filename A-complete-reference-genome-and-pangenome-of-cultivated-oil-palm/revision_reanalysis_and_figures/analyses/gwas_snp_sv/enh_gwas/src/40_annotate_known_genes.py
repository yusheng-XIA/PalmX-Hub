#!${DATA_DIR}/miniconda3/bin/python
"""Map every trait's corrected SNP/SV signals to nearby genes and literature markers."""
from __future__ import annotations

import bisect
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
PREV = Path("${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/04_integrated_SNP_SV_gene_validation_20260721/tables")
LOCUS_GAP, FLANK = 250_000, 100_000
CHROMS = [f"chr{i:02d}B" for i in range(1, 17)]


def known_marker(row):
    seed = str(row.seed_ortholog); desc = str(row.description).lower(); pref = str(row.preferred_name).upper(); pfam = str(row.PFAMs)
    if "XP_008781978" in seed: return "AGL11/STK-like", "A"
    if pref == "WRKY47": return "WRKY47", "B"
    if pref == "PSBS": return "PsbS", "B"
    if pref == "CINV1" or "neutral invertase" in desc: return "neutral invertase", "B"
    if "pectate lyase" in desc or "polygalacturonase" in desc: return "pectin enzyme", "B"
    if "XP_008785886" in seed: return "ARF17-like", "C"
    if "XP_008783395" in seed or ("pot family" in desc and "PTR2" in pfam): return "NPF/PTR transporter", "C"
    if "udp-glycosyltransferase" in desc or "udp glycosyltransferase" in desc: return "UGT family", "C"
    return "", ""


def trait_relevance(marker, category, trait):
    """Conservative relevance flag; physical proximity alone is not trait validation."""
    t = trait.lower()
    rules = {
        "AGL11/STK-like": category == "yield" and any(x in t for x in ["shell", "nut", "kernel", "fruit", "oil"]),
        "WRKY47": category in {"growth", "yield"} or "leaf" in t,
        "PsbS": category == "photosynthesis",
        "neutral invertase": category == "yield" or any(x in t for x in ["fruit", "flesh", "oil"]),
        "pectin enzyme": category == "yield" and any(x in t for x in ["fruit", "flesh", "width", "length", "thickness", "weight"]),
        "ARF17-like": category in {"growth", "yield"},
        "NPF/PTR transporter": category in {"cold_resistant", "growth"},
        "UGT family": category in {"quality", "yield"},
    }
    return "trait_relevant" if marker and rules.get(marker, False) else ("cross_trait_only" if marker else "not_literature_marked")


def load_evidence():
    curated = RUN / "tables/known_gene_literature_reference.tsv"
    if curated.is_file():
        d = pd.read_csv(curated, sep="\t", dtype=str).fillna("")
        return {r.literature_marker: {k: getattr(r, k, "") for k in d.columns} for r in d.itertuples(index=False)}
    frames = []
    for name in ["snp_literature_evidence_marked.tsv", "sv_trait_gene_literature_marked.tsv"]:
        p = PREV / name
        if p.is_file(): frames.append(pd.read_csv(p, sep="\t", dtype=str))
    evidence = {}
    for d in frames:
        for r in d.itertuples(index=False):
            marker = getattr(r, "literature_marker", "") or getattr(r, "reported_gene", "")
            if marker and marker not in evidence:
                evidence[marker] = {k: getattr(r, k, "") for k in ["reported_species", "species", "reported_trait", "evidence_note", "basis", "interpretation", "citation", "doi", "url"]}
    return evidence


def cluster(points):
    result = []
    for chrom in CHROMS:
        sub = sorted([x for x in points if x[1] == chrom], key=lambda x: x[2]); cur = []; anchor = None
        for rec in sub:
            if anchor is None or rec[2] - anchor <= LOCUS_GAP:
                if anchor is None: anchor = rec[2]
                cur.append(rec)
            else:
                result.append(cur); cur = [rec]; anchor = rec[2]
        if cur: result.append(cur)
    return result


def nearby(genes, chrom, start, end):
    g = genes.get(chrom)
    if g is None or g.empty: return g
    hit = g[(g.start <= end + FLANK) & (g.end >= max(1, start - FLANK))]
    if not hit.empty: return hit
    center = (start + end) // 2
    dist = np.where(center < g.start, g.start - center, np.where(center > g.end, center - g.end, 0))
    return g.iloc[[int(np.argmin(dist))]]


def signal_sets(task):
    category, trait = task.category, task.trait
    sdir = RUN / "snp/results" / category / trait; vdir = RUN / "sv/results" / category / trait
    snp_sig = pd.read_csv(sdir / "significant_snps.tsv", sep="\t")
    if len(snp_sig):
        pts = [(r.SNP, r.CHR, int(r.BP), float(r.P), "genomewide_bonferroni") for r in snp_sig.itertuples(index=False)]
    else:
        top = pd.read_csv(sdir / "top1000_snps.tsv", sep="\t").iloc[0]
        pts = [(top.SNP, top.CHR, int(top.BP), float(top.P), "top_signal_exploratory")]
    snp_loci = []
    for i, group in enumerate(cluster(pts), 1):
        lead = min(group, key=lambda x: x[3]); snp_loci.append(("SNP", f"{trait}_SNP_L{i:03d}", group[0][1], min(x[2] for x in group), max(x[2] for x in group), lead[0], lead[2], lead[3], lead[4], len(group)))
    sv_sig = pd.read_csv(vdir / "significant_svs.tsv", sep="\t")
    if len(sv_sig):
        sv_pts = [(r.SV, r.CHR, int(r.BP), float(r.P), "genomewide_bonferroni") for r in sv_sig.itertuples(index=False)]
        sv_rows = []
        for i, group in enumerate(cluster(sv_pts), 1):
            lead = min(group, key=lambda x: x[3])
            sv_rows.append(("SV", f"{trait}_SV_L{i:03d}", group[0][1], min(x[2] for x in group), max(x[2] for x in group),
                            lead[0], lead[2], lead[3], lead[4], len(group)))
    else:
        top = pd.read_csv(vdir / "top1000_svs.tsv", sep="\t").iloc[0]
        sv_rows = [("SV", f"{trait}_SV_L001", top.CHR, int(top.BP), int(top.BP), top.SV, int(top.BP), float(top.P), "top_signal_exploratory", 1)]
    return snp_loci + sv_rows


def main():
    tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
    genes_df = pd.read_csv(RUN / "tables/africa_hap2_gene_catalog.tsv", sep="\t").fillna("")
    genes = {c: g.sort_values("start").reset_index(drop=True) for c, g in genes_df.groupby("chrom")}
    evidence = load_evidence(); loci, candidates = [], []
    for task in tasks.itertuples(index=False):
        for modality, locus_id, chrom, start, end, lead, lead_pos, lead_p, level, nsignals in signal_sets(task):
            loci.append({"category": task.category, "trait": task.trait, "modality": modality, "locus_id": locus_id, "chrom": chrom,
                         "start": start, "end": end, "lead_variant": lead, "lead_pos": lead_pos, "lead_p": lead_p, "signal_level": level, "n_significant_variants": nsignals})
            for gene in nearby(genes, chrom, start, end).itertuples(index=False):
                marker, tier = known_marker(gene)
                distance = 0 if gene.start <= end and gene.end >= start else min(abs(start - gene.end), abs(gene.start - end))
                e = evidence.get(marker, {})
                candidates.append({"category": task.category, "trait": task.trait, "modality": modality, "locus_id": locus_id, "chrom": chrom,
                    "locus_start": start, "locus_end": end, "lead_variant": lead, "lead_pos": lead_pos, "lead_p": lead_p, "signal_level": level,
                    "gene_id": gene.gene_id, "transcript_id": gene.transcript_id, "gene_start": gene.start, "gene_end": gene.end,
                    "relation": "overlap" if distance == 0 else "within_100kb", "distance_bp": distance, "description": gene.description,
                    "preferred_name": gene.preferred_name, "seed_ortholog": gene.seed_ortholog, "GOs": gene.GOs, "KEGG_ko": gene.KEGG_ko,
                    "KEGG_pathway": gene.KEGG_pathway, "PFAMs": gene.PFAMs, "known_gene_marker": marker, "evidence_tier": tier,
                    "biological_relevance": trait_relevance(marker, task.category, task.trait),
                    "reported_species": e.get("reported_species", e.get("species", "")), "reported_trait": e.get("reported_trait", ""),
                    "evidence_note": e.get("evidence_note", e.get("basis", e.get("interpretation", ""))), "citation": e.get("citation", ""), "doi": e.get("doi", ""), "url": e.get("url", "")})
    loci = pd.DataFrame(loci); cand = pd.DataFrame(candidates)
    loci.to_csv(RUN / "tables/all_trait_gwas_loci.tsv", sep="\t", index=False)
    cand.to_csv(RUN / "tables/all_trait_gwas_candidate_genes.tsv", sep="\t", index=False)
    known = cand[cand.known_gene_marker != ""].copy(); known.to_csv(RUN / "tables/literature_supported_known_gene_overlaps.tsv", sep="\t", index=False)
    known[known.biological_relevance == "trait_relevant"].to_csv(RUN / "tables/trait_relevant_known_gene_overlaps.tsv", sep="\t", index=False)
    summaries = []
    for task in tasks.itertuples(index=False):
        x = loci[loci.trait == task.trait]; k = known[known.trait == task.trait]
        kr = k[k.biological_relevance == "trait_relevant"]
        summaries.append({"category": task.category, "trait": task.trait,
            "snp_genomewide_loci": int(((x.modality == "SNP") & (x.signal_level == "genomewide_bonferroni")).sum()),
            "sv_genomewide_loci": int(((x.modality == "SV") & (x.signal_level == "genomewide_bonferroni")).sum()),
            "candidate_gene_rows": int((cand.trait == task.trait).sum()), "known_gene_overlap_rows": len(k),
            "trait_relevant_known_gene_rows": len(kr),
            "known_gene_markers": ";".join(sorted(set(k.known_gene_marker))),
            "trait_relevant_known_gene_markers": ";".join(sorted(set(kr.known_gene_marker)))})
    summary = pd.DataFrame(summaries); summary.to_csv(RUN / "tables/trait_known_gene_overlap_summary.tsv", sep="\t", index=False)

    shown = summary[(summary.snp_genomewide_loci > 0) | (summary.sv_genomewide_loci > 0) | (summary.trait_relevant_known_gene_rows > 0)].copy()
    shown["total"] = shown.snp_genomewide_loci + shown.sv_genomewide_loci
    top = shown.sort_values(["total", "trait_relevant_known_gene_rows"], ascending=False).head(24).reset_index(drop=True)
    y = np.arange(len(top)); fig, (ax, tax) = plt.subplots(1, 2, figsize=(14.2, max(6.8, len(top) * .38)), gridspec_kw={"width_ratios": [1.35, 1]})
    sx = np.log10(top.snp_genomewide_loci.to_numpy() + 1); vx = np.log10(top.sv_genomewide_loci.to_numpy() + 1)
    ax.scatter(sx, y - .13, s=40, color="#365F91", label="SNP windows", zorder=3)
    ax.scatter(vx, y + .13, s=40, color="#7B4F9D", label="SV windows", zorder=3)
    for i, r in top.iterrows():
        if r.snp_genomewide_loci: ax.text(np.log10(r.snp_genomewide_loci + 1) + .035, i - .13, str(r.snp_genomewide_loci), va="center", fontsize=7, color="#365F91")
        if r.sv_genomewide_loci: ax.text(np.log10(r.sv_genomewide_loci + 1) + .035, i + .13, str(r.sv_genomewide_loci), va="center", fontsize=7, color="#7B4F9D")
    ax.set_yticks(y, top.trait); ax.invert_yaxis(); ax.set_xlabel("Bonferroni reporting windows (log10[N+1])")
    ax.set_xticks(np.log10(np.array([1, 2, 11, 101, 1001])), ["0", "1", "10", "100", "1000"])
    ax.set_title("Genome-wide significant windows", loc="left", fontsize=11.5)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="x", color="#E6E6E6", lw=.5); ax.legend(frameon=False, ncol=2, loc="lower right")
    tax.set_ylim(ax.get_ylim()); tax.set_xlim(0, 1); tax.set_xticks([]); tax.set_yticks([]); tax.spines[:].set_visible(False)
    tax.set_title("Trait-relevant literature markers", loc="left", fontsize=11.5)
    for i, r in top.iterrows():
        label = r.trait_relevant_known_gene_markers or "—"
        tax.text(.01, i, textwrap.fill(label.replace(";", "; "), width=48), va="center", fontsize=7.5,
                 color="#8B1A1A" if label != "—" else "#999999")
    fig.suptitle("Zero- and phenotype-tail-excluded GWAS: significant windows and literature-supported nearby candidates\nWindow counts are not LD-independent loci; A/B/C evidence tiers are provided in the companion table.", fontsize=12.5, y=.995)
    fig.tight_layout(); out = RUN / "figures/integrated"; out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "all_trait_known_gene_overlap_overview.png", dpi=420); fig.savefig(out / "all_trait_known_gene_overlap_overview.pdf"); plt.close(fig)
    print(f"traits={len(tasks)} loci={len(loci)} candidate_gene_rows={len(cand)} known_gene_rows={len(known)}")


if __name__ == "__main__": main()
