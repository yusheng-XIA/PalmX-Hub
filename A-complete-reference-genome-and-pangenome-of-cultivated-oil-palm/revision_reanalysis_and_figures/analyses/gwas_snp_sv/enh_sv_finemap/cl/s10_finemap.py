"""Reciprocal conditional analysis and joint SNP+SV Wakefield-ABF fine-mapping for SV-GWAS and SNP-GWAS peaks.

Loci
  SV-anchored : the 115 published Bonferroni SV reporting intervals (8 traits). SV lead = published lead SV;
                SNP lead = top SNP (published SNP model) within +-500 kb of the SV lead.
  SNP-anchored: the published Bonferroni SNP reporting intervals of the same 8 traits plus nut weight (SHELL peak).
                SNP lead = published lead SNP; SV lead = top SV (published SV model) within +-500 kb of the SNP lead.
Models (published EMMAX settings, re-implemented in emmaxpy; REML re-estimated for every conditional fit)
  SNP model: SNP kinship, intercept.   SV model: SV kinship, intercept + SV PC1-5.
Conditional tests
  (a) SNP lead in SNP model + SV lead dosage;  (b) SV lead in SV model + SNP lead dosage   [primary, as published]
  (a2) SNP lead in SV model + SV lead;         (b2) SV lead in SNP model + SNP lead          [common-model check]
  'abolished' = conditional -log10P <= 20% of unconditional AND conditional P above the modality Bonferroni threshold.
Fine-mapping: all SNPs and SVs within +-500 kb of the anchor lead, single-causal-variant Wakefield ABF with equal
prior per variant; prior effect SD = f x SD(trait), f = 0.2 (primary; 0.1 and 0.4 sensitivity); effects and SEs from
the SNP model (primary) or SV model (sensitivity); 95% credible set by cumulative posterior.
Dosages are coded as minor-allele counts among the analysed accessions; missing genotypes mean-imputed (as EMMAX)."""
import sys, time, numpy as np, pandas as pd
sys.path.insert(0, "${CLUSTER_WORK}/enh_gwas/cl")
from config import RUN, BED, FAM, KIN_SNP, SV_TPED, KIN_SV, W as WG
from emmaxpy import GLS, read_bed_rows, read_matrix, read_pheno
from pathlib import Path
O = Path("${CLUSTER_WORK}/enh_sv_finemap"); (O / "out").mkdir(parents=True, exist_ok=True)
SNP_THR, SV_THR, FLANK = 1.7735e-9, 1.3509e-7, 500_000
FRUIT = ["Flesh_thickness_mm", "Nut_length_mm", "Shell_thickness_mm", "Shell_weight_g", "Nut_weight_g"]
t0 = time.time()
ids = [l.split()[1] for l in open(FAM)]; N = len(ids)
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t"); cat_of = dict(zip(tasks.trait, tasks.category))
loci = pd.read_csv(RUN / "tables/all_trait_gwas_loci.tsv", sep="\t")
gw = loci[loci.signal_level == "genomewide_bonferroni"]
SVT = sorted(gw[gw.modality == "SV"].trait.unique())
svl = gw[gw.modality == "SV"].copy(); snpl = gw[(gw.modality == "SNP") & gw.trait.isin(SVT + ["Nut_weight_g"])].copy()
far = pd.read_csv(WG / "out/sv_far_intervals.tsv", sep="\t")
far_pub = set(far.loc[far.far_published, "sv_locus"])
Ks, Kv = read_matrix(KIN_SNP), read_matrix(KIN_SV)
Xs = np.ones((N, 1)); Xv = pd.read_csv(WG / "in/X_sv5pc.tsv", sep="\t", header=None, index_col=0).loc[ids].to_numpy()
rng = pd.read_csv(WG / "in/chr_ranges.tsv", sep="\t").set_index("chr")
POS = {c: np.load(WG / f"in/pos_{c}.npy") for c in rng.index}
# SV genotypes (count of allele '2'; unique ids, first occurrence as in the published per-SV table)
code = {"1": 0, "2": 1}
sid, sch, spos, SG, seen = [], [], [], [], set()
with open(SV_TPED) as fh:
    for line in fh:
        f = line.split()
        if f[1] in seen: continue
        seen.add(f[1])
        a = np.array([code.get(x, -9) for x in f[4:]], dtype=np.int16).reshape(-1, 2)
        d = a.sum(1).astype(np.float32); d[(a < 0).any(1)] = np.nan
        sid.append(f[1]); sch.append(f[0]); spos.append(int(f[3])); SG.append(d)
SG = np.vstack(SG); sid = np.array(sid); sch = np.array(sch); spos = np.array(spos)
sidx = {s: i for i, s in enumerate(sid)}
print("SVs", len(sid), f"{time.time() - t0:.0f}s", flush=True)
meta = pd.read_csv("${ANALYSIS_DIR}/14_pan_genome/06_Minigraph/Pangenie/02_sv_combined/"
                   "步骤四_过滤分类统计/sv_type_stats/sv_qc.per_sv.tsv", sep="\t").drop_duplicates("id").set_index("id")

PH, GL = {}, {}
def pheno(tr):
    if tr not in PH:
        y = read_pheno(RUN / f"phenotypes/{cat_of[tr]}/{tr}.txt", ids); PH[tr] = (y, np.isfinite(y))
    return PH[tr]
def gls(tr, mod):
    if (tr, mod) not in GL:
        y, k = pheno(tr); X, K = (Xs, Ks) if mod == "snp" else (Xv, Kv)
        GL[(tr, mod)] = GLS(y[k], X[k], K[np.ix_(k, k)])
    return GL[(tr, mod)]
def minor(G):
    """orient rows to minor-allele counts; returns oriented copy, MAF, MAC"""
    G = np.array(G, dtype=np.float64); mu = np.nanmean(G, 1); flip = mu > 1
    G[flip] = 2 - G[flip]; nn = (~np.isnan(G)).sum(1); ac = np.nansum(G, 1)
    return G, ac / (2 * nn), ac
def cond(tr, mod, test, cov):
    y, k = pheno(tr); X, K = (Xs, Ks) if mod == "snp" else (Xv, Kv)
    c = np.where(np.isnan(cov), np.nanmean(cov[k]), cov)
    g = GLS(y[k], np.c_[X, c][k], K[np.ix_(k, k)])
    return g.scan(test[k][None, :])[2][0]
def r2(a, b):
    ok = ~np.isnan(a) & ~np.isnan(b)
    if ok.sum() < 3 or a[ok].std() == 0 or b[ok].std() == 0: return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1] ** 2)
def nl(p): return -np.log10(max(float(p), 1e-300))
def abolished(p0, p1, thr): return bool(nl(p1) <= 0.2 * nl(p0) and p1 >= thr)
def classify(snp_ab, sv_ab):
    return {(True, False): "SV_explains_SNP", (False, True): "SNP_explains_SV", (True, True): "mutual", (False, False): "neither"}[(snp_ab, sv_ab)]
def abf(beta, se, sd, f):
    ok = np.isfinite(se) & (se > 0)
    lab = np.full(len(beta), -np.inf); V = se[ok] ** 2; Wp = (f * sd) ** 2; z = beta[ok] / se[ok]
    lab[ok] = 0.5 * np.log(V / (V + Wp)) + 0.5 * z ** 2 * Wp / (V + Wp)
    pp = np.exp(lab - lab.max()); return pp / pp.sum()
def credset(pp, is_sv):
    o = np.argsort(-pp); cs = o[: int(np.searchsorted(np.cumsum(pp[o]), 0.95) + 1)]
    return cs, int(is_sv[cs].sum()), int((~is_sv[cs]).sum())

rows, shell_rows = [], []
def analyse(anchor, tr, locus, chrom, apos, snp_lead_pos=None, sv_lead=None, in_far=False):
    y, k = pheno(tr); n = int(k.sum()); sd = float(np.std(y[k], ddof=1))
    pos = POS[chrom]; a, b = np.searchsorted(pos, [apos - FLANK, apos + FLANK + 1])
    st = int(rng.loc[chrom, "start"])
    gs = read_bed_rows(BED, N, st + a, st + b) if b > a else np.zeros((0, N), np.float32)
    vi = np.nonzero((sch == chrom) & (spos >= apos - FLANK) & (spos <= apos + FLANK))[0]
    if sv_lead is not None and sidx[sv_lead] not in set(vi): vi = np.r_[vi, sidx[sv_lead]]
    if len(vi) == 0:
        rows.append(dict(anchor=anchor, trait=tr, locus_id=locus, chrom=chrom, anchor_pos=apos, n=n, n_snp_window=int(b - a), n_sv_window=0,
                         snp_lead=f"{chrom}:{snp_lead_pos}", snp_lead_pos=snp_lead_pos, class_primary="no_SV_in_window",
                         class_commonmodel_SV="no_SV_in_window", class_commonmodel_SNP="no_SV_in_window"))
        return
    gv = SG[vi]
    Gs, mafs, macs = minor(gs[:, k]); Gv, mafv, macv = minor(gv[:, k])
    fullS = np.array(gs, dtype=np.float64); fullV = np.array(gv, dtype=np.float64)       # full-sample, oriented as Gs/Gv
    fs = np.nanmean(gs[:, k], 1) > 1; fullS[fs] = 2 - fullS[fs]
    fv = np.nanmean(gv[:, k], 1) > 1; fullV[fv] = 2 - fullV[fv]
    res = {}
    for mod in ["snp", "sv"]:
        g = gls(tr, mod)
        bs, ss, ps = g.scan(Gs) if len(Gs) else (np.zeros(0),) * 3
        bv, sv_, pv = g.scan(Gv)
        res[mod] = (np.r_[bs, bv], np.r_[ss, sv_], np.r_[ps, pv])
    ns = len(Gs); is_sv = np.r_[np.zeros(ns, bool), np.ones(len(Gv), bool)]
    names = np.r_[[f"{chrom}:{p}" for p in pos[a:b]], sid[vi]].astype(object)
    vpos = np.r_[pos[a:b], spos[vi]]
    mac_all = np.r_[macs, macv]; maf_all = np.r_[mafs, mafv]
    # leads
    if sv_lead is None:
        j = int(np.argmin(np.where(is_sv, res["sv"][2], np.inf))); sv_lead = names[j]
    jv = int(np.nonzero(names == sv_lead)[0][0])
    if ns == 0:
        js = None
    elif snp_lead_pos is None:
        js = int(np.argmin(res["snp"][2][:ns]))
    else:
        js = int(np.searchsorted(pos[a:b], snp_lead_pos)); assert pos[a + js] == snp_lead_pos
    m = meta.loc[sv_lead] if sv_lead in meta.index else None
    r = dict(anchor=anchor, trait=tr, locus_id=locus, chrom=chrom, anchor_pos=apos, far_from_SNP_interval=in_far, n=n, trait_sd=sd,
             n_snp_window=ns, n_sv_window=int(len(vi)), sv_lead=sv_lead, sv_lead_pos=int(vpos[jv]),
             sv_type=m.svtype if m is not None else "", sv_len=int(m.svlen) if m is not None else np.nan,
             sv_ref_len=int(m.ref_len) if m is not None else np.nan, sv_alt_len=int(m.alt_len) if m is not None else np.nan,
             sv_lead_MAF=maf_all[jv], sv_lead_MAC=int(mac_all[jv]),
             sv_lead_P_SVmodel=res["sv"][2][jv], sv_lead_P_SNPmodel=res["snp"][2][jv],
             sv_lead_beta_SD_SNPmodel=res["snp"][0][jv] / sd, sv_lead_se_SD_SNPmodel=res["snp"][1][jv] / sd)
    if js is not None:
        s_g, v_g = fullS[js], fullV[jv - ns]
        r.update(snp_lead=names[js], snp_lead_pos=int(vpos[js]), snp_lead_MAF=maf_all[js], snp_lead_MAC=int(mac_all[js]),
                 snp_lead_P_SNPmodel=res["snp"][2][js], snp_lead_P_SVmodel=res["sv"][2][js],
                 snp_lead_beta_SD_SNPmodel=res["snp"][0][js] / sd, snp_lead_se_SD_SNPmodel=res["snp"][1][js] / sd,
                 r2_snp_sv=r2(s_g[k], v_g[k]), dist_snp_sv=int(abs(vpos[js] - vpos[jv])))
        r["condA_P_snp_given_sv"] = cond(tr, "snp", s_g, v_g)
        r["condB_P_sv_given_snp"] = cond(tr, "sv", v_g, s_g)
        r["condA2_P_snp_given_sv_SVmodel"] = cond(tr, "sv", s_g, v_g)
        r["condB2_P_sv_given_snp_SNPmodel"] = cond(tr, "snp", v_g, s_g)
        r["snp_retained_frac"] = nl(r["condA_P_snp_given_sv"]) / nl(r["snp_lead_P_SNPmodel"])
        r["sv_retained_frac"] = nl(r["condB_P_sv_given_snp"]) / nl(r["sv_lead_P_SVmodel"])
        sa = abolished(r["snp_lead_P_SNPmodel"], r["condA_P_snp_given_sv"], SNP_THR)
        va = abolished(r["sv_lead_P_SVmodel"], r["condB_P_sv_given_snp"], SV_THR)
        r.update(snp_abolished=sa, sv_abolished=va, class_primary=classify(sa, va))
        sa2 = abolished(r["snp_lead_P_SVmodel"], r["condA2_P_snp_given_sv_SVmodel"], SNP_THR)
        va2 = abolished(r["sv_lead_P_SNPmodel"], r["condB2_P_sv_given_snp_SNPmodel"], SV_THR)
        r["class_commonmodel_SV"] = classify(sa2, va2)
        sa3 = abolished(r["snp_lead_P_SNPmodel"], r["condA_P_snp_given_sv"], SNP_THR)
        r["class_commonmodel_SNP"] = classify(sa3, va2)
        r["lead_snp_stronger_SNPmodel"] = bool(res["snp"][2][js] < res["snp"][2][jv])
        # among all SNP & SV in window, top variant
    else:
        r.update(class_primary="no_SNP_in_window", class_commonmodel_SV="no_SNP_in_window", class_commonmodel_SNP="no_SNP_in_window")
    for mod, f in [("snp", 0.2), ("snp", 0.1), ("snp", 0.4), ("sv", 0.2)]:
        be, se, p = res[mod]; pp = abf(be, se, sd, f); cs, nsv, nsnp = credset(pp, is_sv)
        tag = f"{mod.upper()}model_f{f}"
        top = int(np.argmax(pp))
        r[f"CS95_size_{tag}"] = len(cs); r[f"CS95_nSV_{tag}"] = nsv; r[f"CS95_nSNP_{tag}"] = nsnp
        r[f"svlead_in_CS95_{tag}"] = bool(jv in set(cs.tolist())); r[f"any_SV_in_CS95_{tag}"] = nsv > 0
        r[f"PP_svlead_{tag}"] = pp[jv]; r[f"PP_allSV_{tag}"] = pp[is_sv].sum()
        r[f"PP_snplead_{tag}"] = pp[js] if js is not None else np.nan
        r[f"top_variant_{tag}"] = names[top]; r[f"top_PP_{tag}"] = pp[top]; r[f"top_is_SV_{tag}"] = bool(is_sv[top])
        if f == 0.2 and mod == "snp":
            r["CS95_SVs_SNPmodel_f0.2"] = ";".join(f"{names[i]}({pp[i]:.3f})" for i in cs if is_sv[i])
            best_sv = int(np.argmax(np.where(is_sv, pp, -1))); r["best_SV_by_PP"] = names[best_sv]; r["best_SV_PP"] = pp[best_sv]
            r["n_SNP_P_below_svlead_SNPmodel"] = int((p[:ns] < p[jv]).sum())
    rows.append(r)
    # SHELL region: condition the SNP lead on every SV in the window
    if anchor == "SNP" and chrom == "chr01B" and 3.0e6 <= apos <= 3.5e6 and tr in FRUIT and js is not None:
        s_g = fullS[js]; p0 = res["snp"][2][js]
        for t in range(len(vi)):
            if macv[t] < 1: continue
            pc = cond(tr, "snp", s_g, fullV[t])
            mm = meta.loc[sid[vi[t]]] if sid[vi[t]] in meta.index else None
            shell_rows.append(dict(trait=tr, snp_lead=names[js], snp_lead_P=p0, sv=sid[vi[t]], sv_pos=int(spos[vi[t]]),
                                   sv_type=mm.svtype if mm is not None else "", sv_len=int(mm.svlen) if mm is not None else np.nan,
                                   sv_ref_len=int(mm.ref_len) if mm is not None else np.nan,
                                   sv_MAC=int(macv[t]), sv_P_SNPmodel=res["snp"][2][ns + t], sv_P_SVmodel=res["sv"][2][ns + t],
                                   r2_with_snp_lead=r2(s_g[k], fullV[t][k]), P_snp_lead_given_sv=pc, snp_retained_frac=nl(pc) / nl(p0)))

for x in svl.itertuples(index=False):
    analyse("SV", x.trait, x.locus_id, x.chrom, int(x.lead_pos), sv_lead=x.lead_variant, in_far=x.locus_id in far_pub)
    print(len(rows), x.locus_id, rows[-1]["class_primary"], f"{time.time() - t0:.0f}s", flush=True)
pd.DataFrame(rows).to_csv(O / "out/finemap_SV_anchored.tsv", sep="\t", index=False)
nsv = len(rows)
for x in snpl.itertuples(index=False):
    analyse("SNP", x.trait, x.locus_id, x.chrom, int(x.lead_pos), snp_lead_pos=int(x.lead_pos))
    if len(rows) % 25 == 0: print(len(rows), x.locus_id, rows[-1]["class_primary"], f"{time.time() - t0:.0f}s", flush=True)
D = pd.DataFrame(rows); D.to_csv(O / "out/finemap_all.tsv", sep="\t", index=False)
pd.DataFrame(shell_rows).to_csv(O / "out/shell_snp_lead_conditioned_on_each_SV.tsv", sep="\t", index=False)
print("done", f"{time.time() - t0:.0f}s")
