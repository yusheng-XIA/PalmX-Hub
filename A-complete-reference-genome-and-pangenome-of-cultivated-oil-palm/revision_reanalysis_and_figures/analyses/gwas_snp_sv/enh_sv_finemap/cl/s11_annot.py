"""Annotate lead / best-PP SVs against FL-Hap2 (Africa_hap2) gene models: CDS, UTR (exon minus CDS), intron,
promoter (2 kb upstream of the transcript start, strand-aware), otherwise intergenic with the nearest gene.
SV footprint on FL-Hap2: [pos, pos + ref_len - 1] (DEL, MNV/COMPLEX); insertions [pos, pos + 1].
Function: eggNOG-mapper catalog of FL-Hap2 (seed ortholog, description, preferred name, PFAMs, KEGG)."""
import numpy as np, pandas as pd
from pathlib import Path
O = Path("${CLUSTER_WORK}/enh_sv_finemap")
GFF = "${ANALYSIS_DIR}/20_results/10_database/03_genes/Africa_hap2.gff3"
CAT = ("${ANALYSIS_DIR}/05_GWAS/00_analysis/06_GWAS/"
       "04_integrated_SNP_SV_gene_validation_20260721/tables/africa_hap2_gene_catalog.tsv")
META = ("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/"
        "步骤四_过滤分类统计/sv_type_stats/sv_qc.per_sv.tsv")
g = pd.read_csv(GFF, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 6, 8], names=["chr", "type", "s", "e", "strand", "attr"])
g = g[g.type.isin(["gene", "mRNA", "exon", "CDS"])]
g["gene"] = g.attr.str.extract(r"(evm\.(?:TU|model)\.[^;.]+\.\d+)")[0].str.replace("evm.model.", "evm.TU.", regex=False)
genes = g[g.type == "gene"].set_index("gene")
cat = pd.read_csv(CAT, sep="\t", dtype=str).set_index("gene_id")
meta = pd.read_csv(META, sep="\t").drop_duplicates("id").set_index("id")
D = pd.read_csv(O / "out/finemap_all.tsv", sep="\t")
svs = sorted(set(D.sv_lead.dropna()) | set(D.best_SV_by_PP.dropna()))
def func(gid):
    if gid not in cat.index: return ""
    c = cat.loc[gid]
    parts = [c.preferred_name if c.preferred_name != "-" else "", c.description if c.description != "-" else "",
             f"PFAM:{c.PFAMs}" if c.PFAMs != "-" else "", f"seed:{c.seed_ortholog}"]
    return " | ".join(p for p in parts if p)
out = []
for s in svs:
    m = meta.loc[s]; chrom, pos = m.chrom, int(m.pos)
    e = pos + max(int(m.ref_len), 2) - 1
    sub = g[(g.chr == chrom) & (g.e >= pos - 2000) & (g.s <= e + 2000)]
    hits = {"CDS": set(), "UTR": set(), "intron": set(), "promoter_2kb": set()}
    for gid, q in sub.groupby("gene"):
        gr = genes.loc[gid]
        cds = q[q.type == "CDS"]; ex = q[q.type == "exon"]
        ov = lambda t: ((t.e >= pos) & (t.s <= e)).any()
        if ov(cds): hits["CDS"].add(gid)
        # UTR: exon bases outside CDS overlapped
        if len(ex) and ov(ex):
            cmin, cmax = (cds.s.min(), cds.e.max()) if len(cds) else (np.inf, -np.inf)
            for t in ex.itertuples():
                a, b = max(t.s, pos), min(t.e, e)
                if a <= b and (a < cmin or b > cmax): hits["UTR"].add(gid)
        if gr.s <= e and gr.e >= pos:
            covered = ex[(ex.e >= pos) & (ex.s <= e)]
            span = set()
            # intron if any SV base inside gene lies outside all exons
            lo, hi = max(pos, gr.s), min(e, gr.e)
            in_exon = sum(max(0, min(t.e, hi) - max(t.s, lo) + 1) for t in covered.itertuples())
            if hi - lo + 1 > in_exon: hits["intron"].add(gid)
        ps, pe = (gr.s - 2000, gr.s - 1) if gr.strand == "+" else (gr.e + 1, gr.e + 2000)
        if pe >= pos and ps <= e: hits["promoter_2kb"].add(gid)
    cls = next((k for k in ["CDS", "UTR", "intron", "promoter_2kb"] if hits[k]), "intergenic")
    allg = sorted(set().union(*hits.values()))
    gc = genes[genes.chr == chrom]
    d = np.maximum(0, np.maximum(gc.s - e, pos - gc.e)); j = int(np.argmin(d.to_numpy())); ng = gc.index[j]
    out.append(dict(sv=s, chrom=chrom, pos=pos, end_on_FL_Hap2=e, sv_type=m.svtype, sv_len=int(m.svlen), ref_len=int(m.ref_len),
                    alt_len=int(m.alt_len), feature_class=cls,
                    CDS_genes=",".join(sorted(hits["CDS"])), UTR_genes=",".join(sorted(hits["UTR"])),
                    intron_genes=",".join(sorted(hits["intron"])), promoter_2kb_genes=",".join(sorted(hits["promoter_2kb"])),
                    overlapped_gene_functions=" || ".join(f"{x}: {func(x)}" for x in allg),
                    nearest_gene=ng, nearest_gene_distance_bp=int(d.iloc[j]), nearest_gene_function=func(ng)))
pd.DataFrame(out).to_csv(O / "out/sv_annotation_FL_Hap2.tsv", sep="\t", index=False)
print(pd.DataFrame(out).feature_class.value_counts())
