#!/usr/bin/env python3
"""Independent recomputation of Figure 2c/2d/2e from upstream raw quantification.
Runs on ${COMPUTE_HOST}. Reads only (user dirs read-only); writes to work/trace/Fig2/out/."""
import os, sys, gzip, math
from collections import defaultdict
import numpy as np, pandas as pd

B = "${ANALYSIS_DIR}"
RUN = B + "/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001/outputs"
OUT = "${CLUSTER_WORK}/trace/Fig2/out"
os.makedirs(OUT, exist_ok=True)
STAGES = ["0d","15d","35d","50d","65d","80d","95d","110d","125d","140d","155d","170d","185d","12h","24h","36h","48h","60h","72h"]
def log(*a):
    print(*a, flush=True)

cw = pd.read_csv(RUN + "/stage1_identity/sample_crosswalk_114.tsv", sep="\t", dtype=str)
assert len(cw) == 114
KEYS = cw.integration_key.tolist()

# ---------------------------------------------------------------- metabolome
part = sys.argv[1] if len(sys.argv) > 1 else "all"
if part in ("all", "met"):
    targets = pd.read_csv(RUN + "/stage10_full_coverage_targets_attempt001/full_within_ion_group_representatives.tsv", sep="\t", low_memory=False)
    log("targets", targets.shape, targets.compound_group_id.nunique())
    mem = pd.read_csv(RUN + "/stage26_full_chemical_phenotype_atlas_attempt002/phenotype_axis_compound_memberships.tsv", sep="\t", low_memory=False)
    fig_cg = ["CG06979"]  # placeholder, replaced below from representative table
    rep = pd.read_csv(OUT + "/../Figure2e_representative_metabolites_19stage.tsv", sep="\t")
    fig_cg = rep.compound_group_id.unique().tolist()
    axes = ["P01", "P02", "P03", "P04"]
    mm = mem[mem.axis_id.isin(axes)].copy()
    need = set(fig_cg) | set(mm.compound_group_id)
    t = targets[targets.compound_group_id.isin(need)].copy()
    mats = {}
    for pol in ("NEG", "POS"):
        m = pd.read_csv(RUN + f"/stage2_matched_metabolome_attempt002/{pol.lower()}_matched_114_log2.tsv", sep="\t", index_col=0)
        mats[pol] = m
        log(pol, m.shape)
    rows = []
    for r in t.itertuples():
        v = mats[r.polarity].loc[r.feature_id, KEYS].astype(float).values
        rows.append(v)
    replog = pd.DataFrame(np.vstack(rows), index=t.within_ion_group_id.values, columns=KEYS)
    # z per representative (ddof=1 over finite samples)
    def zrow(v):
        f = np.isfinite(v); n = f.sum()
        if n < 2: return np.full_like(v, np.nan), True
        mu = v[f].mean(); s = np.sqrt(((v[f] - mu) ** 2).sum() / max(n - 1, 1))
        if not np.isfinite(s) or s == 0: return np.full_like(v, np.nan), True
        return (v - mu) / s, False
    Z = {}; BAD = {}
    for wid in replog.index:
        Z[wid], BAD[wid] = zrow(replog.loc[wid].values)
    t["q_rsd"] = t.qc_rsd_after_percent.where(np.isfinite(t.qc_rsd_after_percent), np.inf)
    t["q_kme"] = t.assigned_kME.abs().where(np.isfinite(t.assigned_kME), -np.inf)
    t["nk"] = t.network_keep.astype(str).str.upper().isin(["TRUE", "1", "YES"]).astype(int)
    cons = {}; method = {}
    for cg, g in t.groupby("compound_group_id"):
        g2 = g.sort_values(["qc_detected_fraction", "q_rsd", "n_detected_samples", "q_kme", "nk", "polarity"],
                           ascending=[False, True, False, False, False, True])
        prim = g2.iloc[0].within_ion_group_id
        if len(g) == 1:
            cons[cg] = Z[prim]; method[cg] = "single"
        else:
            a, b = g.within_ion_group_id.tolist()
            ra = pd.Series(replog.loc[a].values); rb = pd.Series(replog.loc[b].values)
            ok = ra.notna() & rb.notna()
            rho = ra[ok].rank().corr(rb[ok].rank())
            if np.isfinite(rho) and rho >= 0.30 and not (BAD[a] or BAD[b]):
                with np.errstate(all="ignore"):
                    cons[cg] = np.nanmean(np.vstack([Z[a], Z[b]]), axis=0)
                method[cg] = "mean2"
            else:
                cons[cg] = Z[prim]; method[cg] = "primary_fallback"
    CZ = pd.DataFrame(cons, index=KEYS).T
    # compare with upstream consensus z matrix rows
    up = pd.read_csv(RUN + "/stage13_compound_abundance_attempt003/compound_consensus_z_matrix.tsv", sep="\t", index_col=0, usecols=None)
    up = up.loc[CZ.index, KEYS]
    d = (CZ - up).abs()
    log("consensus z recompute vs stage13: compounds", len(CZ), "max abs diff", np.nanmax(d.values),
        "NA mismatch", int((CZ.isna() != up.isna()).values.sum()))
    # Fig 2c: 8 metabolites, FL-TN delta of stage means
    recs = []
    for cg in fig_cg:
        name = rep.loc[rep.compound_group_id == cg, "display_name"].iloc[0]
        for s in STAGES:
            vals = {}
            for gt in ("FL", "TN"):
                cols = [f"{gt}|{s}|R{i}" for i in (1, 2, 3)]
                v = CZ.loc[cg, cols].astype(float)
                vals[gt] = (v.mean(), int(v.notna().sum()))
            recs.append(dict(compound_group_id=cg, display_name=name, stage=s, method=method[cg],
                             FL_mean=vals["FL"][0], FL_n=vals["FL"][1], TN_mean=vals["TN"][0], TN_n=vals["TN"][1],
                             delta_FL_minus_TN=vals["FL"][0] - vals["TN"][0]))
    pd.DataFrame(recs).to_csv(OUT + "/fig2c_recomputed.tsv", sep="\t", index=False)
    # Fig 2d: axis scores
    srec = []
    cov = []
    for ax in axes:
        m = mm[mm.axis_id == ax].drop_duplicates("compound_group_id", keep="first")
        w = (m.score_direction.astype(float) * m.confidence_weight.astype(float)).values
        mat = CZ.loc[m.compound_group_id].values.astype(float)
        sc = np.nansum(mat * w[:, None], axis=0) / np.abs(w).sum()
        ev = m.selected_evidence_level.value_counts().to_dict()
        cov.append(dict(axis_id=ax, signed_n=len(m), pos=int((m.score_direction > 0).sum()), neg=int((m.score_direction < 0).sum()),
                        evidence=str(ev), evidence_basis=str(m.evidence_basis.str.contains("MS1").sum())))
        for k, v in zip(KEYS, sc):
            gt, s, r = k.split("|")
            srec.append(dict(axis_id=ax, genotype=gt, stage=s, rep=r, score=v))
    sdf = pd.DataFrame(srec)
    summ = sdf.groupby(["axis_id", "genotype", "stage"]).score.agg(["mean", "std", "count"]).reset_index()
    summ["se"] = summ["std"] / np.sqrt(summ["count"])
    summ.to_csv(OUT + "/fig2d_recomputed.tsv", sep="\t", index=False)
    pd.DataFrame(cov).to_csv(OUT + "/fig2d_axis_membership_check.tsv", sep="\t", index=False)
    upsc = pd.read_csv(RUN + "/stage26_full_chemical_phenotype_atlas_attempt002/molecular_phenotype_scores_by_sample.tsv", sep="\t")
    upsc = upsc[upsc.axis_id.isin(axes)].merge(sdf.assign(integration_key=sdf.genotype + "|" + sdf.stage + "|" + sdf.rep), on=["axis_id", "integration_key"])
    log("axis score per-sample recompute vs stage26 max abs diff", (upsc.molecular_phenotype_score - upsc.score).abs().max(), len(upsc))
    log("metabolome done")

# ---------------------------------------------------------------- RNA
if part in ("all", "rna"):
    counts = pd.read_csv(B + "/22_answer_reviews/00_ms/03_V3/03_figure3/00_minipan/03_rnaseq_mapping/count_matrices/pangraphrna_hisat2_graph.gene_counts.tsv", sep="\t", index_col=0)
    log("counts", counts.shape)
    samples = cw.rna_sample.tolist()
    C = counts[samples].astype(float)
    keep = (C >= 10).sum(axis=1) >= 3
    Ck = C[keep]
    with np.errstate(divide="ignore"):
        lg = np.log(Ck.values)
    lgm = lg.mean(axis=1)
    ok = np.isfinite(lgm)
    sf = []
    for j in range(lg.shape[1]):
        col = lg[:, j]; sel = ok & (Ck.values[:, j] > 0)
        sf.append(float(np.exp(np.median(col[sel] - lgm[sel]))))
    sf = pd.Series(sf, index=samples)
    sfu = pd.read_csv(B + "/22_answer_reviews/00_ms/05_MS/MS_revision_3/Final_20260823_VectorRevision/8月29日/10_Table1_ST8_AI_Cleanup_20260902/ST9_latest_reanalysis_audit/RNA_114_joint_size_factors.tsv", sep="\t").set_index("sample").DESeq2_size_factor
    log("genes kept for SF", int(keep.sum()), "size factor max rel diff vs upstream", float(((sf - sfu.loc[samples]) / sfu.loc[samples]).abs().max()))
    src = pd.read_csv(OUT + "/../Figure2f_RNA_latest114_joint_source.tsv", sep="\t")
    genes = {m: g.split(";") for m, g in src[src.genotype == "FL"].groupby("marker").matched_features.first().items()}
    fa = pd.read_csv(B + "/21_MS/03_result/01_omic/FA_gene_list.tsv", sep="\t", dtype=str)
    fam = {"ACCase": "ACCase", "KASIII": "KASIII", "ENR": "FabI (ENR)", "FATA/B": "FATA/B", "LACS": "LACS", "GPAT": "GPAT", "DGAT": "DGAT", "FAD2-like": "FAD2"}
    chk = []
    for mk, gl in genes.items():
        fl = set(fa.loc[fa.Enzyme == fam.get(mk, "?"), "GeneID"]) if mk in fam else set()
        chk.append(dict(marker=mk, n_genes=len(gl), genes=";".join(gl), FA_gene_list_family_n=len(fl),
                        in_FA_list=len(set(gl) & fl), not_in_FA_list=";".join(sorted(set(gl) - fl))))
    pd.DataFrame(chk).to_csv(OUT + "/fig2e_rna_marker_gene_check.tsv", sep="\t", index=False)
    N = C.div(sf, axis=1)
    recs = []
    for mk, gl in genes.items():
        tot = N.loc[gl].sum(axis=0)
        for s in STAGES:
            mv = {}
            for gt in ("FL", "TN"):
                ss = cw.loc[(cw.genotype == gt) & (cw.stage == s), "rna_sample"].tolist()
                assert len(ss) == 3
                mv[gt] = tot.loc[ss].mean()
            recs.append(dict(marker=mk, stage=s, FL=mv["FL"], TN=mv["TN"], FL_over_TN=mv["FL"] / mv["TN"]))
    pd.DataFrame(recs).to_csv(OUT + "/fig2e_rna_recomputed.tsv", sep="\t", index=False)
    log("RNA done")

# ---------------------------------------------------------------- protein
if part in ("all", "prot"):
    P = B + "/22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724"
    mat = pd.read_csv(P + "/02_current_Astral_114/runs/RUN-PROT-ASTRAL-DIRECTLFQ-20260725-001/outputs/current114_directlfq_protein_abundance_integration_keys.tsv", sep="\t", low_memory=False)
    log("directLFQ matrix", mat.shape)
    raw = pd.read_csv(P + "/02_current_Astral_114/runs/RUN-PROT-ASTRAL-DIRECTLFQ-20260725-001/outputs/current114_final_global1pct_directlfq_input.tsv.diann_precursorsprimary.protein_intensities.tsv", sep="\t", low_memory=False)
    log("raw directLFQ protein_intensities", raw.shape)
    # relabel check: every integration-key column equals exactly one raw run column (on common protein order)
    rawi = raw.set_index("protein").reindex(mat.unified_protein_group)
    mi = mat.set_index("unified_protein_group")
    runmap = {}
    rawv = rawi.fillna(0).values
    for k in KEYS:
        v = mi[k].fillna(0).values
        hits = [rawi.columns[j] for j in range(rawv.shape[1]) if np.allclose(rawv[:, j], v, rtol=1e-9, atol=0)]
        runmap[k] = ";".join(hits)
    pd.Series(runmap).to_csv(OUT + "/protein_key_to_run_map.tsv", sep="\t", header=["raw_run"])
    log("keys with exactly 1 matching raw run", sum(1 for v in runmap.values() if v and ";" not in v))
    log("protein groups total", len(mat))
    # --- independent family assignment (re-implemented)
    fa = pd.read_csv(B + "/21_MS/03_result/01_omic/FA_gene_list.tsv", sep="\t", dtype=str)
    fams = ["ACCase", "KASIII", "FabI (ENR)", "FATA/B", "LACS", "GPAT", "FAD2", "DGAT"]
    fam_genes = {f: set(fa.loc[fa.Enzyme == f, "GeneID"].dropna()) for f in fams}
    ann = pd.read_csv(B + "/20_results/Figure2/07_new_figure/05_omic/1.final_counts/GO_annotation/Africa_hap2/Africa_hap2.emapper.annotations",
                      sep="\t", comment="#", header=None, dtype=str, low_memory=False)
    txt = ann[7].fillna("") + " " + ann[8].fillna("")
    fam_genes["OLE16"] = set(ann.loc[txt.str.contains(r"oleosin 16|\bole16\b", case=False, regex=True), 0])
    fam_genes["LOX (loss risk)"] = set(ann.loc[txt.str.contains(r"plant lipoxygenase|\blox9\b", case=False, regex=True), 0])
    oau = pd.read_csv(B + "/22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-ALLELE-CNS-V4-20260723-001/outputs/stage28_proteome_interface_attempt001/six_genome_protein_to_OAU_metrics.tsv.gz",
                      sep="\t", usecols=["allele_unit_id", "gene_id"], dtype=str)
    g2o = oau.dropna().groupby("gene_id").allele_unit_id.agg(set).to_dict()
    olab = defaultdict(set)
    for f, gs in fam_genes.items():
        for g in gs:
            for o in g2o.get(g, ()):
                olab[o].add(f)
    ouni = {o: next(iter(l)) for o, l in olab.items() if len(l) == 1}
    uni = pd.read_csv(P + "/01_reference_database/FL_TN_unified_exact_sequence_nr_map.tsv", sep="\t", usecols=["unified_protein_id", "source_gene_ids"], dtype=str)
    FAD2_CHR08 = {"UFTN046834", "UFTN046835", "UFTN046836", "UFTN046837"}
    ufam = {}
    for r in uni.itertuples():
        gids = [x.split(":", 1)[1] if ":" in x else x for x in str(r.source_gene_ids).split(";")]
        os_ = set().union(*[g2o.get(g, set()) for g in gids]) if gids else set()
        labs = {ouni[o] for o in os_ if o in ouni}
        if len(labs) == 1:
            f = next(iter(labs))
            if f == "FAD2" and r.unified_protein_id not in FAD2_CHR08:
                continue
            ufam[r.unified_protein_id] = f
    gfam = []
    for grp in mat.unified_protein_group.astype(str):
        labs = {ufam[x] for x in grp.split(";") if x in ufam}
        gfam.append(next(iter(labs)) if len(labs) == 1 else None)
    mat["fam"] = gfam
    mat[mat.fam.notna()][["unified_protein_group", "fam"]].to_csv(OUT + "/protein_group_family_recomputed.tsv", sep="\t", index=False)
    log("assigned groups per family", mat.fam.value_counts().to_dict())
    recs = []
    for f in fams + ["OLE16", "LOX (loss risk)"]:
        sub = mat[mat.fam == f]
        for s in STAGES:
            for gt in ("FL", "TN"):
                reps = []
                for i in (1, 2, 3):
                    v = pd.to_numeric(sub[f"{gt}|{s}|R{i}"], errors="coerce")
                    v = v[np.isfinite(v) & (v > 0)]
                    reps.append(v.sum() if len(v) else np.nan)
                reps = np.array(reps, float)
                recs.append(dict(family=f, stage=s, genotype=gt, groups=len(sub), detected_reps=int(np.isfinite(reps).sum()),
                                 mean_detected=np.nanmean(reps) if np.isfinite(reps).any() else np.nan))
    pd.DataFrame(recs).to_csv(OUT + "/fig2e_protein_recomputed.tsv", sep="\t", index=False)
    log("protein done")
