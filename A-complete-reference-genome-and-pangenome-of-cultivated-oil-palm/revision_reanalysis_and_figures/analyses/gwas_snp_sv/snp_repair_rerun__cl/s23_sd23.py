"""SNP parts of Supplementary Data 23 from the repaired-call-set scan (published model M0), following the published
40_annotate_known_genes.py rules (250-kb anchor clustering; genes overlapping or within 100 kb, else nearest)."""
import numpy as np, pandas as pd, importlib.util
from config import *
spec = importlib.util.spec_from_file_location("ann", RUN / "scripts/40_annotate_known_genes.py"); ann = importlib.util.module_from_spec(spec); spec.loader.exec_module(ann)
O = W / "sd23"; O.mkdir(exist_ok=True)
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
genes_df = pd.read_csv(RUN / "tables/africa_hap2_gene_catalog.tsv", sep="\t").fillna("")
genes = {c: g.sort_values("start").reset_index(drop=True) for c, g in genes_df.groupby("chrom")}
ev = pd.read_csv(RUN / "tables/known_gene_literature_reference.tsv", sep="\t", dtype=str).fillna("")
evidence = {r.literature_marker: r._asdict() for r in ev.itertuples(index=False)}
thr = BONF_SNP
sig_rows, loci, cand = [], [], []
for t in tasks.itertuples(index=False):
    D = W / "out/scan/M0" / t.trait
    hits = pd.concat([pd.read_csv(D / f"{c}.hits.tsv", sep="\t").assign(chrom=c) for c in CHROMS], ignore_index=True)
    sig = hits[hits.p < thr].sort_values("p")
    for r in sig.itertuples(index=False):
        sig_rows.append(dict(category=t.category, trait=t.trait, modality="SNP", variant_id=f"{r.chrom}:{r.pos}", chrom=r.chrom, pos=int(r.pos),
                             p=float(r.p), bonferroni=thr, beta=float(r.beta), se=float(r.se)))
    if len(sig):
        pts = [(f"{r.chrom}:{r.pos}", r.chrom, int(r.pos), float(r.p), "genomewide_bonferroni") for r in sig.itertuples(index=False)]
    else:
        best = None
        for c in CHROMS:
            nlp = np.load(D / f"{c}.npy"); j = int(np.argmax(nlp))
            if best is None or nlp[j] > best[0]: best = (float(nlp[j]), c, int(np.load(W / f"in/pos_{c}.npy")[j]))
        pts = [(f"{best[1]}:{best[2]}", best[1], best[2], 10 ** -best[0], "top_signal_exploratory")]
    for i, g in enumerate(ann.cluster(pts), 1):
        lead = min(g, key=lambda x: x[3]); lid = f"{t.trait}_SNP_L{i:03d}"; start, end = min(x[2] for x in g), max(x[2] for x in g)
        loci.append(dict(category=t.category, trait=t.trait, modality="SNP", locus_id=lid, chrom=g[0][1], start=start, end=end,
                         lead_variant=lead[0], lead_pos=lead[2], lead_p=lead[3], signal_level=lead[4], n_significant_variants=len(g)))
        for gene in ann.nearby(genes, g[0][1], start, end).itertuples(index=False):
            marker, tier = ann.known_marker(gene)
            dist = 0 if gene.start <= end and gene.end >= start else min(abs(start - gene.end), abs(gene.start - end))
            e = evidence.get(marker, {})
            cand.append(dict(category=t.category, trait=t.trait, modality="SNP", locus_id=lid, chrom=g[0][1], locus_start=start, locus_end=end,
                             lead_variant=lead[0], lead_pos=lead[2], lead_p=lead[3], signal_level=lead[4], gene_id=gene.gene_id,
                             transcript_id=gene.transcript_id, gene_start=gene.start, gene_end=gene.end,
                             relation="overlap" if dist == 0 else "within_100kb", distance_bp=dist, description=gene.description,
                             preferred_name=gene.preferred_name, seed_ortholog=gene.seed_ortholog, GOs=gene.GOs, KEGG_ko=gene.KEGG_ko,
                             KEGG_pathway=gene.KEGG_pathway, PFAMs=gene.PFAMs, known_gene_marker=marker, evidence_tier=tier,
                             biological_relevance=ann.trait_relevance(marker, t.category, t.trait),
                             reported_species=e.get("reported_species", ""), reported_trait=e.get("reported_trait", ""),
                             evidence_note=e.get("interpretation", ""), citation=e.get("citation", ""), doi=e.get("doi", ""), url=e.get("url", "")))
pd.DataFrame(sig_rows).to_csv(O / "sig_snp.tsv", sep="\t", index=False)
pd.DataFrame(loci).to_csv(O / "loci_snp.tsv", sep="\t", index=False)
pd.DataFrame(cand).to_csv(O / "cand_snp.tsv", sep="\t", index=False)
print("sig", len(sig_rows), "loci", len(loci), "cand", len(cand))
