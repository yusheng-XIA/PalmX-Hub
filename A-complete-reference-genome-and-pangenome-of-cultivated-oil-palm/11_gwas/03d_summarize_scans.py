"""Per model x trait: lambda_GC over all tests, Bonferroni hits, 250-kb reporting intervals (same anchor rule as the
published 40_annotate_known_genes.py), and reproduction checks of M0 against the published EMMAX output."""
import sys, json, numpy as np, pandas as pd
from scipy.stats import chi2
from config import *
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
cat_of = dict(zip(tasks.trait, tasks.category))
GAP = 250_000

def cluster(df):
    out = []
    for c in CHROMS:
        sub = df[df.chrom == c].sort_values("pos"); cur, anchor = [], None
        for r in sub.itertuples(index=False):
            if anchor is None or r.pos - anchor <= GAP:
                anchor = r.pos if anchor is None else anchor; cur.append(r)
            else:
                out.append(cur); cur, anchor = [r], r.pos
        if cur: out.append(cur)
    return out

summ, loci_all, sig_all = [], [], []
for mdir in sorted((W / "out/scan").iterdir()):
    for tdir in sorted(mdir.iterdir()):
        mod, tr = mdir.name, tdir.name
        if not all((tdir / f"{c}.npy").exists() for c in CHROMS):
            print("incomplete", mod, tr); continue
        nlp = np.concatenate([np.load(tdir / f"{c}.npy") for c in CHROMS]).astype(np.float64)
        assert len(nlp) == N_SNP_TESTS
        lam = float(chi2.isf(10 ** -np.median(nlp), 1) / chi2.ppf(0.5, 1))   # median is order-preserving
        hits = pd.concat([pd.read_csv(tdir / f"{c}.hits.tsv", sep="\t").assign(chrom=c) for c in CHROMS], ignore_index=True)
        sig = hits[hits.p < BONF_SNP].copy(); sig["trait"] = tr; sig["model"] = mod
        sig_all.append(sig)
        cl = cluster(sig)
        for i, g in enumerate(cl, 1):
            lead = min(g, key=lambda r: r.p)
            loci_all.append(dict(model=mod, category=cat_of[tr], trait=tr, modality="SNP", locus_id=f"{tr}_SNP_{mod}_L{i:03d}",
                                 chrom=g[0].chrom, start=min(r.pos for r in g), end=max(r.pos for r in g),
                                 lead_variant=f"{lead.chrom}:{lead.pos}", lead_pos=lead.pos, lead_p=lead.p, n_significant_variants=len(g)))
        rj = json.load(open(tdir / "reml.json"))
        summ.append(dict(model=mod, category=cat_of[tr], trait=tr, n=rj["n"], n_fixed=rj["q"], reml_delta=rj["delta"], pseudo_h2=rj["h2"],
                         lambda_gc=lam, n_sig_snps=len(sig), n_intervals=len(cl), n_chrom=sig.chrom.nunique(),
                         min_p=float(10 ** -nlp.max()), top_snp=(f"{hits.loc[hits.p.idxmin(), 'chrom']}:{hits.loc[hits.p.idxmin(), 'pos']}" if len(hits) else "")))
S = pd.DataFrame(summ); S.to_csv(W / "out/summary_models.tsv", sep="\t", index=False)
pd.DataFrame(loci_all).to_csv(W / "out/snp_loci_models.tsv", sep="\t", index=False)
pd.concat(sig_all).to_csv(W / "out/snp_sig_models.tsv", sep="\t", index=False)
# reproduction: M0 vs published
o = pd.read_csv(W / "out/orig_snp_lambda.tsv", sep="\t")
m0 = S[S.model == "M0"].merge(o, on="trait", suffixes=("", "_pub"))
m0["dlambda"] = m0.lambda_gc - m0.lambda_gc_all
m0["ddelta_rel"] = (m0.reml_delta - m0.reml_delta_pub) / m0.reml_delta_pub
print("M0 vs published: max |d lambda|", m0.dlambda.abs().max(), " n_sig identical in", (m0.n_sig_snps == m0.n_sig).sum(), "/", len(m0),
      " max rel d delta", m0.ddelta_rel.abs().max())
print(m0.loc[m0.n_sig_snps != m0.n_sig, ["trait", "n_sig_snps", "n_sig"]].to_string())
m0.to_csv(W / "out/repro_M0_vs_published.tsv", sep="\t", index=False)
rows = []
for tr in (VALID if "norepro" not in sys.argv else []):
    cat = cat_of[tr]; mx = 0; n = 0
    for c in CHROMS:
        mine = np.load(W / f"out/scan/M0/{tr}/{c}.npy").astype(np.float64)
        pub = pd.read_csv(RUN / f"snp/results/{cat}/{tr}/emmax_{c}.ps", sep="\t", header=None, usecols=[3], dtype=np.float64)[3].to_numpy()
        d = np.abs(mine + np.log10(np.clip(pub, 1e-300, 1))); mx = max(mx, float(d.max())); n += len(d)
    rows.append(dict(trait=tr, n_markers=n, max_abs_dlog10P=mx)); print("full-genome check", rows[-1], flush=True)
if rows: pd.DataFrame(rows).to_csv(W / "out/repro_M0_fullgenome.tsv", sep="\t", index=False)
piv = S.pivot(index="trait", columns="model", values=["lambda_gc", "n_sig_snps", "n_intervals", "n_chrom"])
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 100)
print(piv.to_string())
