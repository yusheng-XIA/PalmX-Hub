"""SV intervals lying > 250 kb from any significant SNP interval of the same trait (published 60_graph_pangenome_value.py
rule), recomputed with M0 (published SNP model) and M1 (+5 PCs) SNP intervals. For M1-distant SV intervals: maximum r2
between the SV lead and every SNP within +-1 Mb (trait-phenotyped accessions, pairwise-complete genotypes) and
conditional EMMAX of the SV lead under the published SV model (SV kinship, intercept, SV PC1-5) plus one SNP dosage:
(a) lead SNP of the nearest significant M1 SNP interval on the same chromosome (if none, the M1 top SNP within +-1 Mb);
(b) the SNP with maximum r2 within +-1 Mb."""
import numpy as np, pandas as pd
from config import *
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno, tped_dosage
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t"); cat_of = dict(zip(tasks.trait, tasks.category))
loci = pd.read_csv(RUN / "tables/all_trait_gwas_loci.tsv", sep="\t")
sv = loci[(loci.modality == "SV") & (loci.signal_level == "genomewide_bonferroni")].copy()
snpl = pd.read_csv(W / "out/snp_loci_models.tsv", sep="\t")
pub_snp = loci[(loci.modality == "SNP") & (loci.signal_level == "genomewide_bonferroni")].assign(model="published")
def gap(a0, a1, b0, b1): return max(0, b0 - a1, a0 - b1)
def nearest(r, S):
    q = S[(S.trait == r.trait) & (S.chrom == r.chrom)]
    if q.empty: return np.inf, None
    d = [(gap(r.start, r.end, x.start, x.end), i) for i, x in q.iterrows()]
    dd, i = min(d); return dd, q.loc[i]
res = {}
for name, S in [("published", pub_snp), ("M0", snpl[snpl.model == "M0"]), ("M1", snpl[snpl.model == "M1"])]:
    dist = [nearest(r, S)[0] for r in sv.itertuples(index=False)]
    sv[f"dist_{name}"] = dist
    res[name] = int((np.array(dist) > 250_000).sum())
print("SV intervals", len(sv), "| >250 kb from a significant SNP interval:", res, flush=True)
far = sv[sv.dist_M0 > 250_000].copy()
# SV genotypes of the leads
need = set(far.lead_variant)
svg = {}
with open(SV_TPED) as fh:
    for line in fh:
        f = line.split(maxsplit=4)
        if f[1] in need:
            svg[f[1]] = tped_dosage(f[4].split())
assert len(svg) == len(need), need - set(svg)
Ksv = read_matrix(KIN_SV); Xsv = pd.read_csv(W / "in/X_sv5pc.tsv", sep="\t", header=None, index_col=0).loc[ids].to_numpy()
rng = pd.read_csv(W / "in/chr_ranges.tsv", sep="\t").set_index("chr")
m1sig = pd.read_csv(W / "out/snp_sig_models.tsv", sep="\t"); m1sig = m1sig[m1sig.model == "M0"]
out = []
for r in far.itertuples(index=False):
    tr = r.trait; y = read_pheno(RUN / f"phenotypes/{cat_of[tr]}/{tr}.txt", ids); keep = np.isfinite(y)
    s = svg[r.lead_variant]
    pos = np.load(W / f"in/pos_{r.chrom}.npy"); a, b = np.searchsorted(pos, [r.lead_pos - 1_000_000, r.lead_pos + 1_000_000 + 1])
    G = read_bed_rows(BED, N, int(rng.loc[r.chrom, "start"]) + a, int(rng.loc[r.chrom, "start"]) + b)[:, keep].astype(np.float64)
    sk = s[keep]
    ok = ~np.isnan(G) & ~np.isnan(sk)[None, :]
    n = ok.sum(1); gx = np.where(ok, G, 0); sx = np.where(ok, sk[None, :], 0)
    mg = gx.sum(1) / n; ms = sx.sum(1) / n
    cov = (gx * sx).sum(1) / n - mg * ms
    vg = (gx ** 2).sum(1) / n - mg ** 2; vs = (sx ** 2).sum(1) / n - ms ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where((vg > 0) & (vs > 0), cov ** 2 / (vg * vs), 0)
    ok_sv = ~np.isnan(sk); acs = np.nansum(sk); mac_sv = int(min(acs, 2 * ok_sv.sum() - acs))
    a250, b250 = np.searchsorted(pos, [r.start - 250_000, r.end + 250_000 + 1])
    base = GLS(y[keep], Xsv[keep], Ksv[np.ix_(keep, keep)]).scan(s[keep][None, :])[2][0]
    common = dict(trait=tr, sv_locus=r.locus_id, chrom=r.chrom, sv_start=r.start, sv_end=r.end, sv_lead=r.lead_variant,
                  sv_lead_p_published=r.lead_p, sv_lead_p_recomputed=base, n=int(keep.sum()), sv_lead_insample_mac=mac_sv,
                  dist_to_published_snp_interval=r.dist_published, dist_to_M0_snp_interval=r.dist_M0,
                  n_snps_250kb=int(b250 - a250), n_snps_1Mb=int(b - a))
    if b - a == 0:
        out.append(dict(**common, max_r2=np.nan)); print(out[-1], flush=True); continue
    j = int(np.nanargmax(r2)); p_idx = a + j
    nlp1 = np.load(W / f"out/scan/M0/{tr}/{r.chrom}.npy")[a:b]
    # covariate SNPs
    dnear, near = nearest(r, snpl[snpl.model == "M0"])
    if near is not None:
        cpos = int(near.lead_pos); ctype = "nearest significant SNP interval lead"
    else:
        k = int(np.argmax(nlp1)); cpos = int(pos[a + k]); ctype = "top SNP within 1 Mb (no significant SNP interval on chromosome)"
    ci = int(np.searchsorted(pos, cpos)); assert pos[ci] == cpos
    def cond(gpos_idx):
        g = read_bed_rows(BED, N, int(rng.loc[r.chrom, "start"]) + gpos_idx, int(rng.loc[r.chrom, "start"]) + gpos_idx + 1)[0].astype(np.float64)
        g = np.where(np.isnan(g), np.nanmean(g[keep]), g)
        X = np.c_[Xsv, g][keep]
        gl = GLS(y[keep], X, Ksv[np.ix_(keep, keep)])
        return gl.scan(s[keep][None, :])[2][0]
    out.append(dict(**common, max_r2=float(r2[j]), max_r2_snp=f"{r.chrom}:{pos[p_idx]}", max_r2_snp_M1_p=float(10 ** -nlp1[j]),
                    cond_a_snp=f"{r.chrom}:{cpos}", cond_a_type=ctype, cond_a_snp_r2=float(r2[ci - a]) if a <= ci < b else np.nan,
                    cond_a_p=cond(ci), cond_b_snp=f"{r.chrom}:{pos[p_idx]}", cond_b_p=cond(p_idx),
                    top_M1_snp_1Mb=f"{r.chrom}:{pos[a + int(np.argmax(nlp1))]}", top_M1_snp_1Mb_p=float(10 ** -nlp1.max()),
                    top_M1_snp_1Mb_r2=float(r2[int(np.argmax(nlp1))]), cond_c_p=cond(a + int(np.argmax(nlp1)))))
    print(out[-1], flush=True)
O = pd.DataFrame(out)
O["snp_covered"] = O.n_snps_250kb > 0
O["far_M0"] = O.dist_to_M0_snp_interval > 250_000
O["cond_a_sig"] = O.cond_a_p < BONF_SV; O["cond_b_sig"] = O.cond_b_p < BONF_SV; O["cond_c_sig"] = O.cond_c_p < BONF_SV
O.to_csv(W / "out/sv_far_intervals_M0.tsv", sep="\t", index=False)
sv.to_csv(W / "out/sv_intervals_distances_M0.tsv", sep="\t", index=False)
for lab in ["far_M0"]:
    print(lab, "median max r2", O[O[lab]].max_r2.median())
    x = O[O[lab]]; c = x[x.snp_covered]
    print(f"{lab}: {len(x)} SV intervals; no SNP tested within 250 kb: {(~x.snp_covered).sum()}; SNP-covered: {len(c)};"
          f" max r2 < 0.2: {(c.max_r2 < 0.2).sum()}; max r2 < 0.5: {(c.max_r2 < 0.5).sum()}; still < SV Bonferroni conditional on nearest SNP lead: {c.cond_a_sig.sum()};"
          f" conditional on max-r2 SNP: {c.cond_b_sig.sum()}; conditional on top SNP within 1 Mb: {c.cond_c_sig.sum()}; all three: {(c.cond_a_sig & c.cond_b_sig & c.cond_c_sig).sum()}; a&b: {(c.cond_a_sig & c.cond_b_sig).sum()}; r2<0.2 and both: {((c.max_r2 < 0.2) & c.cond_a_sig & c.cond_b_sig).sum()}")
    print(x.groupby("trait").agg(n=("sv_locus", "size"), uncovered=("snp_covered", lambda v: (~v).sum()), r2lt02=("max_r2", lambda v: (v < 0.2).sum()),
                                 ca=("cond_a_sig", "sum"), cb=("cond_b_sig", "sum"), cc=("cond_c_sig", "sum")).to_string())
