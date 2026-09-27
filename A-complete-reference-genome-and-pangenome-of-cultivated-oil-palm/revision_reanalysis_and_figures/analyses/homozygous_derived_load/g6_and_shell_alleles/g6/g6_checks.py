"""G6 robustness checks on the per-site table (g6_sites_full.tsv.gz) and genotypes (g6_gt.npz).
Per accession: derived alleles, homozygous-derived genotypes, heterozygous genotypes and called sites per class.
Contrast: commercial source groups (SEA-A + SA-EG, n = 105) vs AFR + IDB + SEA-B (n = 108).
Effect = median(COMM) / median(NONC) - 1; chromosome jackknife s.e.; accession bootstrap 95% CI."""
import sys, numpy as np, pandas as pd
from scipy import stats
H = "${WORK_DIR}/fix/headline_pop/g6"
S = pd.read_csv(f"{H}/g6_sites_full.tsv.gz", sep="\t", low_memory=False)
z = np.load(f"{H}/g6_gt.npz"); GT = z["GT"].astype(np.int8); samples = list(z["samples"])
grp = {l.split()[0]: l.split()[1].split(",") for l in open(f"{H}/../../enh_sn_pop/cl/groups308.txt")}
arch = np.array([grp[s][0] for s in samples]); k4 = np.array([grp[s][1] for s in samples])
C = np.isin(arch, ["SEAA", "SAEG"]); N = np.isin(arch, ["AFR", "IDB", "SEAB"])
print("sites", len(S), "samples", len(samples), "COMM", C.sum(), "NONC", N.sum())
dupmask = S.duplicated(["Chrom", "Pos"], keep=False).to_numpy()
S = S[~dupmask].reset_index(drop=True); GT = GT[~dupmask]
print("dropped duplicated-position records", int(dupmask.sum()))
S["callrate"] = (GT >= 0).mean(1)
S["NigerianOK"] = S.Nigerian_aligned.astype(str) == "True"
se_map = {"missense_variant": "missense", "synonymous_variant": "synonymous", "stop_gained": "stop_gained"}
S["se_cat"] = S.SnpEff_effect.map(lambda e: next((v for k, v in se_map.items() if isinstance(e, str) and e.split("&")[0] == k), "other"))
mz = np.where((S.MZ4_h1 == S.MZ4_h2) & S.MZ4_h1.isin(["REF", "ALT"]), S.MZ4_h1, "NA")
S["ole"] = mz
S["cds_frac"] = (S.CDS_idx + 1) / S.CDS_len
S["last_exon"] = S.Exon_idx == S.N_exons - 1
rng = np.random.default_rng(7)


def derived_dose(anc, sel):
    """dosage of the derived allele for selected sites given ancestral state REF/ALT"""
    g = GT[sel].astype(float); g[g < 0] = np.nan
    a = anc[sel]
    d = np.where((a == "REF")[:, None], g, 2 - g)
    return d


def per_sample(d, chrom):
    out = {}
    for c in list(np.unique(chrom)) + [None]:
        m = np.ones(len(chrom), bool) if c is None else chrom != c
        x = d[m]
        out[c] = (np.nansum(x, 0), np.nansum(x == 2, 0), np.nansum(x == 1, 0), np.sum(~np.isnan(x), 0))
    return out


def effects(anc_col, site_filter, cats=("missense", "synonymous", "stop_gained"), cat_col="Cat", tag="", boot=0, keep=False):
    anc = S[anc_col].to_numpy()
    base = (S.callrate >= 0.8) & np.isin(anc, ["REF", "ALT"]) & site_filter
    res = {}; ps = {}
    for cat in cats:
        sel = (base & (S[cat_col] == cat)).to_numpy()
        d = derived_dose(anc, sel)
        ps[cat] = (per_sample(d, S.Chrom.to_numpy()[sel]), int(sel.sum()))
    chroms = sorted(set(S.Chrom))
    def eff(stat_fn, c=None, idxC=None, idxN=None):
        v = stat_fn(c if c is None or all(c in ps[k][0] for k in cats) else None)
        a = v[C] if idxC is None else v[C][idxC]; b = v[N] if idxN is None else v[N][idxN]
        return np.median(a) / np.median(b) - 1
    rows = []
    for cat in cats:
        P, n = ps[cat]
        for lab, j in (("derived_alleles", 0), ("hom_derived", 1)):
            fn = lambda c, j=j, P=P: P[c][j]
            full = eff(fn); jk = np.array([eff(fn, c) for c in chroms]); k = len(jk)
            jse = np.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum())
            bt = []
            for _ in range(boot):
                bt.append(eff(fn, None, rng.integers(0, C.sum(), C.sum()), rng.integers(0, N.sum(), N.sum())))
            lo, hi = (np.percentile(bt, [2.5, 97.5]) if boot else (np.nan, np.nan))
            p = stats.mannwhitneyu(fn(None)[C], fn(None)[N]).pvalue
            rows.append(dict(check=tag, cat=cat, n_sites=n, measure=lab, COMM_median=np.median(fn(None)[C]), NONC_median=np.median(fn(None)[N]),
                             rel_diff=full, jk_se=jse, boot_lo=lo, boot_hi=hi, P_MWU=p))
    # ratios: missense/synonymous and stop/synonymous for hom derived and derived alleles
    for num in [c for c in cats if c != "synonymous"]:
        for lab, j in (("hom_ratio_to_syn", 1), ("derived_ratio_to_syn", 0)):
            fn = lambda c, j=j, num=num: ps[num][0][c][j] / ps["synonymous"][0][c][j]
            full = eff(fn); jk = np.array([eff(fn, c) for c in chroms]); k = len(jk)
            jse = np.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum())
            bt = [eff(fn, None, rng.integers(0, C.sum(), C.sum()), rng.integers(0, N.sum(), N.sum())) for _ in range(boot)]
            lo, hi = (np.percentile(bt, [2.5, 97.5]) if boot else (np.nan, np.nan))
            p = stats.mannwhitneyu(fn(None)[C], fn(None)[N]).pvalue
            rows.append(dict(check=tag, cat=num, n_sites=ps[num][1], measure=lab, COMM_median=np.median(fn(None)[C]), NONC_median=np.median(fn(None)[N]),
                             rel_diff=full, jk_se=jse, boot_lo=lo, boot_hi=hi, P_MWU=p))
    if keep:
        return pd.DataFrame(rows), ps
    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 300)
    allsites = pd.Series(True, index=S.index)
    out = []
    # 0. primary (date palm, own classes) with accession bootstrap
    prim, ps = effects("Phoenix", allsites, tag="primary_datepalm", boot=2000, keep=True); out.append(prim)
    # 1. SnpEff concordance and SnpEff classes
    ph = (S.Phoenix.isin(["REF", "ALT"])) & (S.callrate >= 0.8)
    ct = pd.crosstab(S.Cat[ph], S.se_cat[ph]); print("own class x SnpEff (Phoenix-polarized, callrate>=0.8)\n", ct)
    seav = S.SnpEff_effect.notna() & (S.SnpEff_effect.astype(str) != "NA")
    print("SnpEff available (5 chromosomes):", int((ph & seav).sum()), "polarized sites; concordance by own class:",
          {c: round(float(((S.Cat == c) & (S.se_cat == c) & ph & seav).sum() / ((S.Cat == c) & ph & seav).sum()), 4) for c in ("missense", "synonymous", "stop_gained")})
    out.append(effects("Phoenix", seav, tag="own_classes_5chrom_with_snpeff"))
    out.append(effects("Phoenix", seav, cat_col="se_cat", tag="snpeff_classes_5chrom"))
    agree = (S.Cat == S.se_cat) | ~seav
    out.append(effects("Phoenix", agree, tag="own_and_snpeff_agree"))
    # 2. E. oleifera polarization; both outgroups agree
    out.append(effects("ole", allsites, tag="Eoleifera_polarized"))
    S["both"] = np.where(S.Phoenix == S.ole, S.Phoenix, "NA")
    out.append(effects("both", allsites, tag="datepalm_and_Eoleifera_agree"))
    # 3. reference bias
    out.append(effects("Phoenix", S.NigerianOK, tag="Nigerian_aligned_sites"))
    out.append(effects("Phoenix", S.Phoenix == "REF", tag="derived_is_ALT"))
    out.append(effects("Phoenix", S.Phoenix == "ALT", tag="derived_is_REF"))
    # 4. LoF filters
    lof_ok = ~((S.Cat == "stop_gained") & ((S.cds_frac > 0.95) | S.last_exon))
    out.append(effects("Phoenix", lof_ok, tag="stop_not_last_exon_not_last5pct"))
    out.append(effects("Phoenix", lof_ok & agree, tag="stop_filtered_and_snpeff_agree"))
    R = pd.concat(out); R.to_csv(f"{H}/g6_checks.tsv", sep="\t", index=False)
    print(R[["check", "cat", "n_sites", "measure", "COMM_median", "NONC_median", "rel_diff", "jk_se", "boot_lo", "boot_hi", "P_MWU"]].round(4).to_string())
    # stop-gained site audit
    st = S[(S.Cat == "stop_gained") & S.Phoenix.isin(["REF", "ALT"]) & (S.callrate >= 0.8)]
    print("stop-gained polarized sites", len(st), "last exon", st.last_exon.mean().round(3), "in last 5% of CDS", (st.cds_frac > 0.95).mean().round(3),
          "SnpEff stop_gained", (st.se_cat == "stop_gained").mean().round(3))
    st.to_csv(f"{H}/stop_gained_audit.tsv", sep="\t", index=False)
    # reference-distance covariate: per-sample ALT fraction at synonymous sites
    anc = S.Phoenix.to_numpy(); sel = ((S.callrate >= 0.8) & (S.Cat == "synonymous")).to_numpy()
    g = GT[sel].astype(float); g[g < 0] = np.nan; altfrac = np.nanmean(g, 0) / 2
    P = ps
    df = pd.DataFrame(dict(sample=samples, arch=arch, k4=k4, COMM=C.astype(int), NONC=N.astype(int), altfrac=altfrac))
    for cat in ("missense", "synonymous", "stop_gained"):
        d, h, he, n = P[cat][0][None]
        df[f"{cat}_der"] = d; df[f"{cat}_hom"] = h; df[f"{cat}_het"] = he; df[f"{cat}_n"] = n
    df["het_coding"] = (df.missense_het + df.synonymous_het) / (df.missense_n + df.synonymous_n)
    df.to_csv(f"{H}/per_sample_load_final.tsv", sep="\t", index=False)
    import statsmodels.formula.api as smf
    s = df[(df.COMM == 1) | (df.NONC == 1)].copy()
    s["mh"] = s.missense_hom / s.synonymous_hom; s["sh"] = s.stop_gained_hom / s.synonymous_hom
    for y in ["missense_hom", "mh", "sh"]:
        m = smf.ols(f"{y} ~ COMM + altfrac", data=s).fit()
        print(y, "COMM %.4g P %.2g | altfrac P %.2g" % (m.params.COMM, m.pvalues.COMM, m.pvalues.altfrac))
    print("altfrac COMM %.4f NONC %.4f" % (s[s.COMM == 1].altfrac.median(), s[s.COMM == 0].altfrac.median()))
