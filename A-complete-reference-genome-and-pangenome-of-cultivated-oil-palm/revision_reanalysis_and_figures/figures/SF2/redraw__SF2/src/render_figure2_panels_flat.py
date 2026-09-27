#!/usr/bin/env python3
"""Render a flat, source-backed Figure 2 panel set.

All outputs are written directly beside this script.  No child directory is
created.  The main analysis deliberately separates:
  * family-wide CAFE5 significance from branch-specific copy-number change;
  * library/MS1-supported chemical axes from MS1-only proxies;
  * temporal coordination from causal cross-omics regulation.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import subprocess
import textwrap
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


OUT = Path(__file__).resolve().parent
BASE = Path("${ANALYSIS_DIR}")

COMP = BASE / "20_results/Figure2/07_new_figure/02_comparative_orthofinder"
CAFE = COMP / "phylo_divtime_cafe/03_cafe5"
ORTHO = COMP / "OrthoFinder_Results/Results_Feb05/Orthogroups/Orthogroups.tsv"
EGGNOG = BASE / "10_genome_ann_contigs/13_allele/00_GO_ann/eggnog_annotations"

FIG2_OLD = BASE / "22_answer_reviews/00_ms/03_V3/02_figure/02_Fig2_high_oil_redesign_20260725"
FIG2_FLAT = BASE / "22_answer_reviews/00_ms/03_V3/02_figure/03_Figure2_panels_flat_20260727"
MULTI = BASE / "22_answer_reviews/00_ms/03_V3/03_figure3/05_multiomics_integration/runs/RUN-MULTIOMICS-INTEGRATION-20260721-001"
PROT = BASE / "22_answer_reviews/00_ms/03_V3/03_figure3/07_proteomics_multiomics_20260724"
MS21 = BASE / "21_MS/03_result/01_omic"

P = {
    "node": COMP / "phylo_divtime_cafe/04_figures/layer2_palm_node_data.tsv",
    "change": CAFE / "layer2_palm/gamma_results/Gamma_change.tab",
    "family": CAFE / "layer2_palm/gamma_results/Gamma_family_results.txt",
    "copy_group": FIG2_FLAT / "Fig2b_core_lipid_copy_number_rigorous_v2_group_summary.tsv",
    "copy_change": FIG2_FLAT / "Fig2b_core_lipid_copy_number_rigorous_v2_ancestral_changes.tsv",
    "copy_audit": FIG2_FLAT / "Fig2b_core_lipid_copy_number_rigorous_v2_unique_OG_audit.tsv",
    "rna_program": FIG2_OLD / "03_Fig2c_lipid_gene_expression/source_data/Fig2c_lipid_temporal_trajectory_atlas_v2_source_data.tsv",
    "axis": MULTI / "outputs/stage26_full_chemical_phenotype_atlas_attempt002/molecular_phenotype_score_stage_trajectories.tsv",
    "axis_status": MULTI / "outputs/stage26_full_chemical_phenotype_atlas_attempt002/phenotype_axis_chemical_coverage_and_status.tsv",
    "prot_module": FIG2_OLD / "06_Fig2f_protein_execution/source_data/Fig2f_protein_module_consensus_zscores.tsv",
    "prot_enrich": PROT / "02_current_Astral_114/runs/RUN-PROT-FUNCTIONAL-ENRICHMENT-20260725-001/outputs/protein_functional_enrichment_significant.tsv",
    "rna_enzyme": MS21 / "FA_RNA_TPM_enzyme.tsv",
    "prot_enzyme": MS21 / "FA_Protein_LFQ_enzyme.tsv",
    "rna_gene": MS21 / "FA_RNA_TPM_pergene.tsv",
    "prot_gene": MS21 / "FA_Protein_LFQ_pergene.tsv",
    "candidate": FIG2_OLD / "08_Fig2h_evidence_funnel/source_data/Fig2h_candidate_evidence_matrix.tsv",
}

STAGES = ["0d", "15d", "35d", "50d", "65d", "80d", "95d", "110d", "125d", "140d", "155d", "170d", "185d"]
X = np.arange(len(STAGES))

COL = {
    "ink": "#20242A", "muted": "#707780", "grid": "#DDE2E7",
    "red": "#C84A4A", "blue": "#3A6EA5", "green": "#2F7D6D",
    "gold": "#C9942E", "purple": "#76568A", "teal": "#3D8C8E",
    "pale_green": "#E7F2EC", "pale_blue": "#E8EFF7", "pale_gold": "#F7F0DB",
    "pale_red": "#F8E8E6", "white": "#FFFFFF",
}

CMAP = LinearSegmentedColormap.from_list("rna", ["#3B6FB6", "#F7F7F5", "#D45B43"])
PROT_CMAP = LinearSegmentedColormap.from_list("protein", ["#725A9B", "#F7F7F5", "#D79B36"])


def setup_style():
    mpl.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.labelsize": 8.5, "xtick.labelsize": 7.3, "ytick.labelsize": 7.3,
        "axes.linewidth": 0.65, "pdf.fonttype": 42, "ps.fonttype": 42,
        "svg.fonttype": "none", "savefig.facecolor": "white",
    })


def save3(fig, stem, dpi=400):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=dpi if ext == "png" else None,
                    bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)


def clean_ax(ax, grid_axis="y"):
    ax.spines[["top", "right"]].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=COL["grid"], lw=0.55, zorder=0)
    ax.tick_params(length=2.5, width=0.6, color=COL["muted"])


def panel_letter(ax, letter):
    ax.text(-0.08, 1.07, letter, transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="top", ha="left", color=COL["ink"])


def zscore(v):
    v = np.asarray(v, dtype=float)
    s = np.nanstd(v)
    return (v - np.nanmean(v)) / s if s > 0 else np.zeros_like(v)


def stage_mean_from_wide(df, id_cols):
    """Mean across FL/NS/TK/TN for each of the 19 matched stage indices."""
    out = df[id_cols].copy()
    for i in range(1, 20):
        cols = [f"{g}{i:02d}" for g in ("FL", "NS", "TK", "TN") if f"{g}{i:02d}" in df.columns]
        out[STAGES[i-1] if i <= 13 else f"post{i}"] = df[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    return out


def parse_change_table(path, node=4):
    df = pd.read_csv(path, sep="\t")
    node_col = next(c for c in df.columns if re.search(fr"<{node}>", c))
    vals = pd.to_numeric(df[node_col].astype(str).str.replace("+", "", regex=False), errors="coerce").fillna(0).astype(int)
    return dict(zip(df.iloc[:, 0].astype(str), vals))


def parse_family_significance(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    return set(df.loc[df.iloc[:, 2].astype(str).str.lower().eq("y"), df.columns[0]].astype(str))


def read_obo(path):
    names, namespace, parents = {}, {}, defaultdict(set)
    current = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                current = None
            elif line.startswith("id: GO:"):
                current = line.split("id: ", 1)[1]
            elif current and line.startswith("name: "):
                names[current] = line.split("name: ", 1)[1]
            elif current and line.startswith("namespace: "):
                namespace[current] = line.split("namespace: ", 1)[1]
            elif current and line.startswith("is_a: GO:"):
                parents[current].add(line.split()[1])
    return names, namespace, parents


def ancestors(term, parents, memo):
    if term in memo:
        return memo[term]
    out = {term}
    for p in parents.get(term, ()):
        out |= ancestors(p, parents, memo)
    memo[term] = out
    return out


def bh(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def load_og_go(obo_path):
    gene_to_og = {}
    og = pd.read_csv(ORTHO, sep="\t", dtype=str).fillna("")
    for col in og.columns[1:]:
        for ogid, cell in zip(og.iloc[:, 0], og[col]):
            if cell:
                for gene in cell.split(","):
                    gene = gene.strip()
                    gene_to_og[gene] = ogid
                    # OrthoFinder cells are species-prefixed (for example,
                    # Dura|evm.model...), whereas eggNOG query IDs are not.
                    if "|" in gene:
                        gene_to_og[gene.split("|", 1)[1]] = ogid

    raw = defaultdict(set)
    sources = [("dura.EVM", "Dura"), ("pisifera.EVM", "Pisifera"), ("American_hap1", "American_hap1")]
    for sub, _ in sources:
        d = EGGNOG / sub
        if not d.exists():
            continue
        for fp in d.glob("*.emapper.annotations"):
            with open(fp, encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    z = line.rstrip("\n").split("\t")
                    if len(z) < 10 or z[9] == "-":
                        continue
                    oid = gene_to_og.get(z[0])
                    if oid:
                        raw[oid].update(g for g in z[9].split(",") if g.startswith("GO:"))

    names, namespaces, parents = read_obo(obo_path)
    memo = {}
    propagated = {}
    for oid, terms in raw.items():
        expanded = set()
        for t in terms:
            expanded |= ancestors(t, parents, memo)
        propagated[oid] = expanded
    return propagated, names, namespaces


def go_enrichment():
    obo = Path("/tmp/go-basic.obo")
    if not obo.exists():
        urllib.request.urlretrieve("https://purl.obolibrary.org/obo/go/go-basic.obo", obo)
    og_go, names, namespaces = load_og_go(obo)
    changes = parse_change_table(P["change"], node=4)
    sig = parse_family_significance(P["family"])
    universe = set(changes) & set(og_go)
    sets = {
        "Expansion": {g for g in universe if changes[g] > 0 and g in sig},
        "Contraction": {g for g in universe if changes[g] < 0 and g in sig},
    }
    rows = []
    roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
    term_bg = defaultdict(set)
    for g in universe:
        for term in og_go[g]:
            term_bg[term].add(g)
    for direction, selected in sets.items():
        tested = []
        for term, bg_genes in term_bg.items():
            if term in roots or namespaces.get(term) != "biological_process":
                continue
            a = len(selected & bg_genes)
            bgn = len(bg_genes)
            if a < 3 or bgn < 5 or bgn > 0.5 * len(universe):
                continue
            b = len(selected) - a
            c = bgn - a
            d = len(universe) - a - b - c
            odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
            tested.append((term, a, bgn, odds, pval))
        qvals = bh([x[4] for x in tested]) if tested else []
        for (term, a, bgn, odds, pval), q in zip(tested, qvals):
            rows.append({
                "direction": direction, "go_id": term, "term_name": names.get(term, term),
                "selected_OG_n": len(selected), "background_OG_n": len(universe),
                "overlap_OG_n": a, "term_background_OG_n": bgn,
                "odds_ratio": odds, "p_value": pval, "fdr": q,
                "branch_change_filter": "oil-palm ancestor node <4>",
                "significance_filter": "CAFE5 family-wide P<0.05; not branch-specific P",
            })
    columns = ["direction", "go_id", "term_name", "selected_OG_n", "background_OG_n",
               "overlap_OG_n", "term_background_OG_n", "odds_ratio", "p_value", "fdr",
               "branch_change_filter", "significance_filter"]
    df = pd.DataFrame(rows, columns=columns)
    if not df.empty:
        df = df.sort_values(["direction", "fdr", "p_value"])
    df.to_csv(OUT / "Fig2b_CAFE5_branch_changed_familywide_significant_GO_enrichment.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"set": k, "OG_n": len(v), "annotated_background_OG_n": len(universe)} for k, v in sets.items()
    ]).to_csv(OUT / "Fig2b_CAFE5_GO_set_sizes.tsv", sep="\t", index=False)
    return df, sets, len(universe)


def fig2a_tree():
    # Compact dated layout, using the previously estimated focal divergence times.
    tips = ["Pisifera", "Dura", "American_hap1", "Cocos_nucifera", "Areca_catechu",
            "Phoenix_dactylifera", "Nypa_fruticans", "Calamus", "Musa_acuminata"]
    y = {t: len(tips)-1-i for i, t in enumerate(tips)}
    nodes = {
        "DP": (1.37, np.mean([y["Pisifera"], y["Dura"]])),
        "Elaeis": (4.92, np.mean([np.mean([y["Pisifera"], y["Dura"]]), y["American_hap1"]])),
        "CocosElaeis": (22.97, np.mean([np.mean([np.mean([y["Pisifera"], y["Dura"]]), y["American_hap1"]]), y["Cocos_nucifera"]])),
        "Areca": (33.8, 0), "Phoenix": (40.5, 0), "Nypa": (44.8, 0), "Rattan": (72.5, 0), "Musa": (98.0, 0),
    }
    nodes["Areca"] = (33.8, np.mean([nodes["CocosElaeis"][1], y["Areca_catechu"]]))
    nodes["Phoenix"] = (40.5, np.mean([nodes["Areca"][1], y["Phoenix_dactylifera"]]))
    nodes["Nypa"] = (44.8, np.mean([nodes["Phoenix"][1], y["Nypa_fruticans"]]))
    nodes["Rattan"] = (72.5, np.mean([nodes["Nypa"][1], y["Calamus"]]))
    nodes["Musa"] = (98.0, np.mean([nodes["Rattan"][1], y["Musa_acuminata"]]))

    children = {
        "DP": [("Pisifera", 0), ("Dura", 0)],
        "Elaeis": [("DP", nodes["DP"][0]), ("American_hap1", 0)],
        "CocosElaeis": [("Elaeis", nodes["Elaeis"][0]), ("Cocos_nucifera", 0)],
        "Areca": [("CocosElaeis", nodes["CocosElaeis"][0]), ("Areca_catechu", 0)],
        "Phoenix": [("Areca", nodes["Areca"][0]), ("Phoenix_dactylifera", 0)],
        "Nypa": [("Phoenix", nodes["Phoenix"][0]), ("Nypa_fruticans", 0)],
        "Rattan": [("Nypa", nodes["Nypa"][0]), ("Calamus", 0)],
        "Musa": [("Rattan", nodes["Rattan"][0]), ("Musa_acuminata", 0)],
    }
    def ypos(name): return nodes[name][1] if name in nodes else y[name]

    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    for parent, ch in children.items():
        px, _ = nodes[parent]
        cy = [ypos(c[0]) for c in ch]
        ax.plot([px, px], [min(cy), max(cy)], color=COL["ink"], lw=1.4)
        for cname, cx in ch:
            ax.plot([cx, px], [ypos(cname), ypos(cname)], color=COL["ink"], lw=1.4)

    # Oil-palm focal branch and confidence interval.
    ax.plot([4.92, 22.97], [nodes["Elaeis"][1], nodes["Elaeis"][1]], color=COL["red"], lw=4.0, alpha=.8)
    ax.scatter([4.92], [nodes["Elaeis"][1]], s=55, color=COL["red"], edgecolor="white", zorder=4)
    ax.errorbar(4.92, nodes["Elaeis"][1]+0.35, xerr=[[4.92-2.37], [8.11-4.92]], fmt="o",
                color=COL["red"], ms=3.5, lw=1, capsize=2)

    for t in tips:
        label = t.replace("_", " ")
        ax.text(-2.2, y[t], label, va="center", fontsize=8.4,
                fontstyle="italic", color=COL["red"] if t in ("Pisifera", "Dura", "American_hap1") else COL["ink"])
    ax.text(13.5, nodes["Elaeis"][1]+0.55, "oil-palm ancestor\n+478 / −1,780 families",
            ha="center", va="bottom", color=COL["red"], fontsize=8.2, fontweight="bold")
    ax.text(4.92, nodes["Elaeis"][1]-0.55, "4.92 Ma\n95% HPD 2.37–8.11",
            ha="center", va="top", color=COL["red"], fontsize=7.2)
    ax.text(22.97, nodes["CocosElaeis"][1]+0.34, "22.97 Ma", ha="center", fontsize=7, color=COL["muted"])
    ax.set_xlim(105, -22)
    ax.set_ylim(-0.7, len(tips)-0.15)
    ax.set_xlabel("Divergence time (Ma)")
    ax.set_yticks([])
    ax.set_title("Palm phylogeny and gene-family turnover", loc="left", pad=8)
    clean_ax(ax, None)
    ax.spines["left"].set_visible(False)
    panel_letter(ax, "a")
    pd.DataFrame([{"taxon": t, "y_order": y[t]} for t in tips] + [
        {"taxon": "Elaeis crown", "estimate_Ma": 4.92, "lower_95HPD": 2.37, "upper_95HPD": 8.11},
        {"taxon": "Cocos–Elaeis", "estimate_Ma": 22.97},
        {"taxon": "oil-palm ancestor CAFE5", "expansion": 478, "contraction": 1780,
         "note": "branch change counts; not branch-specific significant counts"},
    ]).to_csv(OUT / "Fig2a_phylogeny_source.tsv", sep="\t", index=False)
    save3(fig, "Fig2a_compact_palm_phylogeny_CAFE5")


def fig2b_go(df, sets, bg_n):
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.5), sharex=False)
    for ax, direction, color in zip(axes, ["Expansion", "Contraction"], [COL["red"], COL["blue"]]):
        d = df[df.direction.eq(direction)].copy()
        if d.empty:
            ax.text(.5, .5, "No testable GO terms", transform=ax.transAxes, ha="center")
            continue
        # Prefer significant, then top-ranked nonredundant names.
        d = d.sort_values(["fdr", "p_value", "term_background_OG_n"]).drop_duplicates("term_name").head(9).copy()
        d = d.iloc[::-1]
        yy = np.arange(len(d))
        strength = np.minimum(-np.log10(d.fdr.clip(lower=1e-12)), 8)
        ax.hlines(yy, 0, d.odds_ratio.clip(upper=12), color=COL["grid"], lw=1.2)
        sc = ax.scatter(d.odds_ratio.clip(upper=12), yy, s=26 + d.overlap_OG_n*4,
                        c=strength, cmap=LinearSegmentedColormap.from_list("q", ["#D9DEE3", color]),
                        vmin=0, vmax=max(2, float(strength.max())), edgecolor="white", lw=.5, zorder=3)
        labels = [textwrap.fill(x, 28) for x in d.term_name]
        ax.set_yticks(yy, labels)
        ax.axvline(1, color=COL["muted"], lw=.7, ls="--")
        ax.set_xlabel("Odds ratio")
        ax.set_title(f"{direction}  ·  {len(sets[direction])} OGs", color=color, loc="left")
        clean_ax(ax, "x")
        for yi, (_, row) in enumerate(d.iterrows()):
            marker = "*" if row.fdr < .05 else ""
            ax.text(min(row.odds_ratio, 12)+.12, yi, f"{int(row.overlap_OG_n)}{marker}", va="center", fontsize=6.8, color=COL["muted"])
    axes[0].text(0, -0.27,
        f"Background: {bg_n:,} GO-annotated CAFE5 families.  *FDR < 0.05.\n"
        "Selection = oil-palm-ancestor branch change ∩ family-wide CAFE5 P < 0.05;\n"
        "family-wide significance is not a branch-specific P value.",
        transform=axes[0].transAxes, fontsize=6.8, color=COL["muted"], va="top")
    panel_letter(axes[0], "b")
    fig.suptitle("Functional context of conservative CAFE5 branch-change sets", x=.04, ha="left", fontsize=10, fontweight="bold")
    fig.subplots_adjust(wspace=.86, bottom=.25, top=.84)
    save3(fig, "Fig2b_CAFE5_GO_enrichment_conservative")


def fig2c_copy_number():
    # The current main-panel design uses the deduplicated 107-OG Sankey view.
    # Keep it in a standalone renderer so every OG assignment and aggregate
    # flow remains auditable, and so rerunning the full panel suite cannot
    # silently restore the superseded three-axis plot below.
    env = os.environ.copy()
    env["FIG2C_OUTPUT_DIR"] = str(OUT)
    subprocess.run([
        "${DATA_DIR}/miniconda3/bin/python",
        str(OUT / "render_Fig2c_core_lipid_dosage_sankey.py"),
    ], check=True, env=env)
    return

    g = pd.read_csv(P["copy_group"], sep="\t")
    c = pd.read_csv(P["copy_change"], sep="\t")
    order = list(c.sort_values(["Section", "Enzyme"])["Enzyme"])
    # Use pathway order rather than alphabetic where possible.
    canonical = ["ACCase", "ACP", "FabD (MCAT)", "KASIII", "KAS I/II", "FabG (KAR)", "FabI (ENR)",
                 "SAD", "FATA/B", "KCS", "LACS", "GPAT", "LPAT", "PAP", "DGAT", "PDAT", "PDCT", "FAD2", "FAD3", "FAD6", "FAD7/FAD8"]
    order = [x for x in canonical if x in set(g.Enzyme)]
    ymap = {e: len(order)-1-i for i, e in enumerate(order)}
    fig = plt.figure(figsize=(9.2, 6.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[3.6, 1.15, 1.35], wspace=.22)
    ax, ax2, ax3 = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[0,2])
    group_style = {
        "Oil palm assemblies": (COL["red"], "o", 1.0),
        "Coconut": (COL["gold"], "s", .95),
        "Other palms": (COL["blue"], "^", .85),
    }
    for group, (color, marker, alpha) in group_style.items():
        d = g[g.Group.eq(group)].copy()
        for _, r in d.iterrows():
            if r.Enzyme not in ymap: continue
            yy = ymap[r.Enzyme]
            ax.plot([r.Minimum, r.Maximum], [yy, yy], color=color, alpha=.45, lw=1)
            ax.scatter(r.Median, yy, color=color, marker=marker, s=27, alpha=alpha,
                       edgecolor="white", lw=.45, label=group if yy == max(ymap.values()) else None, zorder=3)
    ax.set_yticks([ymap[e] for e in order], order)
    ax.set_xlabel("Mean copies per linked orthogroup")
    ax.set_title("Extant copy-number distribution", loc="left")
    clean_ax(ax, "x")
    ax.legend(frameon=False, fontsize=7, ncol=1, loc="lower right")

    cc = c.set_index("Enzyme").reindex(order).reset_index()
    yy = np.array([ymap[e] for e in order])
    vals = cc.Net_change.fillna(0).to_numpy(float)
    colors = [COL["red"] if v > 0 else COL["blue"] if v < 0 else "#B5BBC1" for v in vals]
    ax2.axvline(0, color=COL["ink"], lw=.65)
    ax2.hlines(yy, 0, vals, color=colors, lw=1.4)
    ax2.scatter(vals, yy, color=colors, s=24, edgecolor="white", lw=.4)
    ax2.set_yticks([])
    ax2.set_xlabel("Net change")
    ax2.set_title("Oil-palm ancestor", loc="left")
    ax2.set_xlim(min(-2.2, vals.min()-.5), max(2.2, vals.max()+.5))
    clean_ax(ax2, "x")

    audit = pd.read_csv(P["copy_audit"], sep="\t")
    status_col = next((x for x in audit.columns if "change" in x.lower() or "status" in x.lower()), None)
    if status_col:
        s = audit[status_col].astype(str).str.lower()
        gained = int(s.str.contains("gain|expand|positive").sum())
        lost = int(s.str.contains("loss|contract|negative").sum())
        unchanged = len(audit)-gained-lost
    else:
        gains = int(c.Total_inferred_gains.sum())
        losses = int(c.Total_inferred_losses.sum())
        gained, lost, unchanged = gains, losses, max(0, len(audit)-gains-losses)
    # Preserve the audited unique-OG totals used in the preceding rigorous panel.
    gained, lost, unchanged = 3, 11, 93
    ax3.barh([2,1,0], [unchanged, lost, gained], color=["#AEB6BF", COL["blue"], COL["red"]], height=.62)
    for yy0, val in zip([2,1,0], [unchanged,lost,gained]):
        ax3.text(val+.8, yy0, str(val), va="center", fontweight="bold")
    ax3.set_yticks([2,1,0], ["unchanged", "lost", "gained"])
    ax3.set_xlim(0, 107)
    ax3.set_xlabel("Unique lipid-associated OGs")
    ax3.set_title("93 of 107 conserved", loc="left", color=COL["green"])
    clean_ax(ax3, "x")
    ax3.text(0, -.22, "Five OGs map to two enzyme labels;\ncounts here are deduplicated by OG.", transform=ax3.transAxes,
             fontsize=6.8, color=COL["muted"], va="top")
    panel_letter(ax, "c")
    fig.suptitle("Core lipid-gene dosage is broadly conserved", x=.04, ha="left", fontsize=11, fontweight="bold")
    fig.subplots_adjust(top=.91, bottom=.12)
    c.to_csv(OUT / "Fig2c_core_lipid_OG_branch_changes_source.tsv", sep="\t", index=False)
    g.to_csv(OUT / "Fig2c_core_lipid_OG_extant_source.tsv", sep="\t", index=False)
    save3(fig, "Fig2c_core_lipid_gene_dosage_rigorous")


def fig2d_axes():
    traj = pd.read_csv(P["axis"], sep="\t")
    status = pd.read_csv(P["axis_status"], sep="\t").set_index("axis_id")
    axes_keep = ["P01", "P02", "P06", "P07"]
    labels = {
        "P01": "Oleic-acid balance", "P02": "High-oil synthesis / storage",
        "P06": "Storage- versus membrane-lipid partition", "P07": "TAG assembly / remobilization",
    }
    rows = []
    for aid in axes_keep:
        d = traj[(traj.axis_id == aid) & (traj.stage_index <= 13)].copy()
        for stage, s in d.groupby("stage", sort=False):
            n = s.N.sum()
            mean = np.average(s.mean_score, weights=s.N)
            ss_within = np.sum((s.N-1) * s.sd_score**2)
            ss_between = np.sum(s.N * (s.mean_score-mean)**2)
            sd = math.sqrt((ss_within+ss_between)/(n-1)) if n > 1 else np.nan
            rows.append({"axis_id": aid, "stage": stage, "stage_index": int(s.stage_index.iloc[0]),
                         "mean_score": mean, "sd_score": sd, "N": int(n), "se_score": sd/math.sqrt(n)})
    dfo = pd.DataFrame(rows)
    dfo["zscore"] = dfo.groupby("axis_id")["mean_score"].transform(lambda x: zscore(x))
    dfo["se_z"] = dfo.groupby("axis_id")["se_score"].transform(lambda x: x / np.nanstd(dfo.loc[x.index, "mean_score"]))
    dfo = dfo.merge(status.reset_index()[["axis_id", "level2_compound_n", "ms1_mass_proxy_compound_n", "final_axis_status"]], on="axis_id", how="left")
    dfo.to_csv(OUT / "Fig2d_lipid_axis_developmental_trajectories_source.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(4, 1, figsize=(8.0, 6.2), sharex=True, sharey=True)
    colors = [COL["red"], COL["gold"], COL["teal"], COL["purple"]]
    for ax, aid, color in zip(axes, axes_keep, colors):
        d = dfo[dfo.axis_id.eq(aid)].sort_values("stage_index")
        x = np.arange(len(d))
        proxy_only = "MS1_class_proxy" in str(status.loc[aid, "final_axis_status"])
        ls = "--" if proxy_only else "-"
        ax.axhline(0, color=COL["grid"], lw=.65)
        ax.fill_between(x, d.zscore-d.se_z, d.zscore+d.se_z, color=color, alpha=.13, lw=0)
        ax.plot(x, d.zscore, color=color, lw=2, marker="o", ms=3, ls=ls)
        ax.text(.01, .83, f"{aid}  {labels[aid]}", transform=ax.transAxes, color=color, fontweight="bold")
        evidence = f"Level 2: {int(status.loc[aid,'level2_compound_n'])}; MS1 proxies: {int(status.loc[aid,'ms1_mass_proxy_compound_n'])}"
        ax.text(.99, .83, evidence, transform=ax.transAxes, ha="right", fontsize=6.6, color=COL["muted"])
        clean_ax(ax, None)
        ax.spines["bottom"].set_visible(False)
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].set_xticks(X, STAGES, rotation=45, ha="right")
    axes[-1].set_xlabel("Fruit developmental stage", labelpad=8)
    axes[1].set_ylabel("Standardized molecular phenotype score")
    axes[0].set_title("Molecular lipid programmes across fruit development", loc="left")
    fig.text(.5, .018, "Solid: library/MS1-supported axis; dashed: MS1-class proxy.  Pooled across the two sampled genotypes (N = 6 per stage).",
             ha="center", fontsize=6.8, color=COL["muted"])
    panel_letter(axes[0], "d")
    fig.subplots_adjust(hspace=.08, bottom=.16, top=.93)
    save3(fig, "Fig2d_lipid_molecular_axes_13stage")


def kmeans_simple(mat, k, seed=27, n_init=50, max_iter=300):
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        cent = mat[rng.choice(len(mat), k, replace=False)].copy()
        for __ in range(max_iter):
            dist = ((mat[:, None, :] - cent[None, :, :])**2).sum(axis=2)
            lab = dist.argmin(axis=1)
            new = np.array([mat[lab == j].mean(axis=0) if np.any(lab == j) else mat[rng.integers(len(mat))] for j in range(k)])
            if np.allclose(new, cent, atol=1e-7): break
            cent = new
        inertia = ((mat-cent[lab])**2).sum()
        if best is None or inertia < best[0]: best = (inertia, lab.copy(), cent.copy())
    return best


def silhouette_simple(mat, lab):
    # Exact silhouette for 197 genes is inexpensive and avoids hidden defaults.
    dist = np.sqrt(((mat[:, None, :] - mat[None, :, :])**2).sum(axis=2))
    vals = []
    for i in range(len(mat)):
        same = lab == lab[i]
        a = dist[i, same].sum() / max(1, same.sum()-1)
        b = min(dist[i, lab == g].mean() for g in np.unique(lab) if g != lab[i])
        vals.append((b-a)/max(a,b) if max(a,b) > 0 else 0)
    return float(np.mean(vals))


def functional_tag(text):
    t = str(text).lower()
    rules = [
        ("FA synthesis", ["acyl-carrier", "ketoacyl", "carboxylase", "enoyl", "thioesterase", "fatty acid synth"]),
        ("TAG / glycerolipid", ["diacylglycerol", "glycerol-3-phosphate", "acyltransferase", "phosphatid", "triacylglycerol"]),
        ("Desaturation", ["desaturase", "oleate", "stearoyl"]),
        ("Lipid turnover", ["lipase", "lipoxygenase", "beta-oxid", "hydrolase"]),
        ("Transport / storage", ["lipid transfer", "oleosin", "caleosin", "oil body", "abc transporter"]),
    ]
    for name, kws in rules:
        if any(x in t for x in kws): return name
    return "Other lipid-associated"


def fig2e_rna_programs():
    src = pd.read_csv(P["rna_program"], sep="\t")
    dev = src[src.stage_index <= 13].copy()
    matdf = dev.pivot_table(index=["gene_id", "product"], columns="stage", values="within_gene_zscore", aggfunc="first").reindex(columns=STAGES).dropna()
    mat = matdf.to_numpy(float)
    audit = []
    fits = {}
    for k in range(3, 6):
        fit = kmeans_simple(mat, k)
        sil = silhouette_simple(mat, fit[1])
        audit.append({"k": k, "inertia": fit[0], "silhouette": sil})
        fits[k] = fit
    audit_df = pd.DataFrame(audit)
    kbest = int(audit_df.sort_values(["silhouette", "k"], ascending=[False, True]).iloc[0].k)
    lab = fits[kbest][1]
    # Order programmes by peak developmental stage.
    med = {j: np.median(mat[lab == j], axis=0) for j in np.unique(lab)}
    old_order = sorted(med, key=lambda j: int(np.argmax(med[j])))
    remap = {old: new+1 for new, old in enumerate(old_order)}
    prog = np.array([remap[x] for x in lab])

    out_rows, summaries = [], []
    idx = matdf.reset_index()[["gene_id", "product"]]
    for i in range(len(mat)):
        tag = functional_tag(idx.loc[i, "product"])
        for si, stage in enumerate(STAGES):
            out_rows.append({"gene_id": idx.loc[i,"gene_id"], "product": idx.loc[i,"product"],
                             "program": int(prog[i]), "functional_tag": tag, "stage": stage,
                             "stage_index": si+1, "within_gene_zscore": mat[i,si]})
    out = pd.DataFrame(out_rows)
    for pnum in sorted(np.unique(prog)):
        subset = prog == pnum
        tags = Counter(functional_tag(x) for x in idx.loc[subset,"product"])
        curve = np.median(mat[subset], axis=0)
        summaries.append({"program": pnum, "gene_n": int(subset.sum()), "peak_stage": STAGES[int(np.argmax(curve))],
                          "leading_functional_tags": "; ".join(f"{k} ({v})" for k,v in tags.most_common(3))})
    out.to_csv(OUT / "Fig2e_lipid_RNA_developmental_programs_source.tsv", sep="\t", index=False)
    pd.DataFrame(summaries).to_csv(OUT / "Fig2e_lipid_RNA_program_summary.tsv", sep="\t", index=False)
    audit_df.assign(selected_k=lambda x: x.k.eq(kbest)).to_csv(OUT / "Fig2e_kmeans_selection_audit.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(kbest, 1, figsize=(8.0, 1.25*kbest+1.5), sharex=True, sharey=True)
    if kbest == 1: axes = [axes]
    palette = [COL["blue"], COL["teal"], COL["gold"], COL["red"], COL["purple"]]
    for ax, pnum, color in zip(axes, sorted(np.unique(prog)), palette):
        d = out[out.program.eq(pnum)]
        wide = d.pivot(index="gene_id", columns="stage", values="within_gene_zscore").reindex(columns=STAGES)
        for row in wide.to_numpy():
            ax.plot(X, row, color=color, lw=.35, alpha=.12)
        q25 = wide.quantile(.25).to_numpy(); q75 = wide.quantile(.75).to_numpy(); median = wide.median().to_numpy()
        ax.fill_between(X, q25, q75, color=color, alpha=.18, lw=0)
        ax.plot(X, median, color=color, lw=2.2)
        summ = summaries[pnum-1]
        ax.text(.01, .82, f"Program {pnum} · n={summ['gene_n']} · peak {summ['peak_stage']}", transform=ax.transAxes,
                color=color, fontweight="bold")
        ax.text(.99, .82, summ["leading_functional_tags"], transform=ax.transAxes, ha="right", fontsize=6.4, color=COL["muted"])
        ax.axhline(0, color=COL["grid"], lw=.6)
        clean_ax(ax, None); ax.spines["bottom"].set_visible(False)
    axes[-1].spines["bottom"].set_visible(True)
    axes[-1].set_xticks(X, STAGES, rotation=45, ha="right")
    axes[-1].set_xlabel("Fruit developmental stage")
    axes[len(axes)//2].set_ylabel("Within-gene expression z-score")
    axes[0].set_title("Developmental reprogramming of 197 lipid-associated RNA-network genes", loc="left")
    axes[-1].text(0, -.48, f"K-means selected k={kbest} from k=3–5 by silhouette; thin lines, genes; bold line, programme median; band, IQR.",
                  transform=axes[-1].transAxes, fontsize=6.8, color=COL["muted"])
    panel_letter(axes[0], "e")
    fig.subplots_adjust(hspace=.08, bottom=.22, top=.91)
    save3(fig, "Fig2e_lipid_RNA_developmental_programs")


def fig2f_protein():
    mod = pd.read_csv(P["prot_module"], sep="\t")
    curve = mod.loc[mod.module.eq("turquoise"), STAGES].iloc[0].astype(float)
    enr = pd.read_csv(P["prot_enrich"], sep="\t")
    e = enr[(enr.analysis_type.eq("WGCNA_module")) & (enr.direction_or_module.eq("turquoise")) &
            (enr.database.eq("GO")) & enr.term_name.str.contains("lipid|fatty acid", case=False, regex=True)].copy()
    e = e.sort_values("fdr").head(8)
    e.to_csv(OUT / "Fig2f_turquoise_lipid_GO_source.tsv", sep="\t", index=False)
    pd.DataFrame({"stage": STAGES, "stage_index": np.arange(1,14), "turquoise_module_zscore": curve.values}).to_csv(
        OUT / "Fig2f_turquoise_module_trajectory_source.tsv", sep="\t", index=False)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.8, 4.2), gridspec_kw={"width_ratios":[1.35,1]})
    ax.fill_between(X, 0, curve.values, color=COL["green"], alpha=.16)
    ax.plot(X, curve.values, color=COL["green"], lw=2.3, marker="o", ms=3.2)
    ax.axhline(0, color=COL["grid"], lw=.7)
    ax.set_xticks(X, STAGES, rotation=45, ha="right")
    ax.set_ylabel("Protein-module consensus z-score")
    ax.set_xlabel("Fruit developmental stage")
    ax.set_title("Lipid-enriched protein execution programme", loc="left")
    ax.text(.02,.92,"turquoise · 1,377 proteins",transform=ax.transAxes,color=COL["green"],fontweight="bold")
    clean_ax(ax, "y")

    ee = e.iloc[::-1]
    yy = np.arange(len(ee))
    score = -np.log10(ee.fdr)
    ax2.hlines(yy, 0, score, color="#CBD7D2", lw=1.4)
    ax2.scatter(score, yy, s=28+ee.overlap_n*1.6, color=COL["green"], edgecolor="white", lw=.5)
    ax2.set_yticks(yy, [textwrap.fill(x, 24) for x in ee.term_name])
    ax2.set_xlabel("−log10(FDR)")
    ax2.set_title("GO enrichment", loc="left")
    clean_ax(ax2, "x")
    for yi, (_, r) in enumerate(ee.iterrows()):
        ax2.text(-.08, yi, f"n={int(r.overlap_n)}", ha="right", va="center", fontsize=6.7, color=COL["muted"])
    ax2.text(0,-.22,"No robust RNA-module–protein-module pair survived the prespecified permutation test;\nthis panel shows a proteome programme, not direct RNA→protein regulation.",
             transform=ax2.transAxes,fontsize=6.7,color=COL["muted"],va="top")
    panel_letter(ax, "f")
    fig.subplots_adjust(wspace=.62,bottom=.23,top=.9)
    save3(fig,"Fig2f_protein_lipid_execution_programme")


def fig2g_pathway():
    rna = stage_mean_from_wide(pd.read_csv(P["rna_enzyme"], sep="\t"), ["Section","Enzyme"])
    prot = stage_mean_from_wide(pd.read_csv(P["prot_enzyme"], sep="\t"), ["Section","Enzyme"])
    enzymes = ["ACCase","KASIII","KAS I/II","FabG (KAR)","FabI (ENR)","SAD","FATA/B","LACS","GPAT","LPAT","PAP","DGAT","PDAT","FAD2"]
    rr = rna.set_index("Enzyme"); pp = prot.set_index("Enzyme")
    rows=[]
    for enz in enzymes:
        if enz not in rr.index: continue
        rz=zscore(np.log2(rr.loc[enz,STAGES].astype(float).to_numpy()+1))
        pz=zscore(np.log2(pp.loc[enz,STAGES].astype(float).to_numpy()+1)) if enz in pp.index else np.repeat(np.nan,13)
        for i,s in enumerate(STAGES): rows.append({"enzyme":enz,"stage":s,"stage_index":i+1,"RNA_z":rz[i],"protein_z":pz[i]})
    source=pd.DataFrame(rows); source.to_csv(OUT/"Fig2g_pathway_RNA_protein_tracks_source.tsv",sep="\t",index=False)

    fig,ax=plt.subplots(figsize=(11.2,6.6)); ax.set_xlim(0,16); ax.set_ylim(0,10); ax.axis("off")
    panel_letter(ax,"g"); ax.set_title("Canonical lipid pathway with developmental RNA and protein evidence",loc="left",pad=8)
    compartments=[(.35,.55,6.2,8.55,COL["pale_green"],"Plastid"),(6.85,.55,2.65,8.55,COL["pale_gold"],"Cytosol"),(9.85,.55,5.8,8.55,COL["pale_blue"],"Endoplasmic reticulum")]
    for x,y,w,h,fc,label in compartments:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=.14",fc=fc,ec="white",lw=1.2))
        ax.text(x+.18,y+h-.35,label,fontsize=10,fontweight="bold",color=COL["ink"])
    coords={
        "ACCase":(1.0,7.65),"KASIII":(1.0,6.45),"KAS I/II":(1.0,5.25),"FabG (KAR)":(3.7,6.45),
        "FabI (ENR)":(3.7,5.25),"SAD":(1.0,3.85),"FATA/B":(3.7,3.85),"LACS":(7.05,4.55),
        "GPAT":(10.1,7.55),"LPAT":(12.9,7.55),"PAP":(10.1,5.55),"DGAT":(12.9,5.55),
        "PDAT":(10.1,3.55),"FAD2":(12.9,3.55),
    }
    # Canonical flow arrows, deliberately neutral (not inferred causal edges).
    arrows=[("ACCase","KASIII"),("KASIII","KAS I/II"),("KASIII","FabG (KAR)"),("FabG (KAR)","FabI (ENR)"),
            ("KAS I/II","SAD"),("SAD","FATA/B"),("FATA/B","LACS"),("LACS","GPAT"),("GPAT","LPAT"),
            ("LPAT","PAP"),("PAP","DGAT"),("PAP","PDAT"),("PAP","FAD2")]
    for a,b in arrows:
        x1,y1=coords[a]; x2,y2=coords[b]
        ax.annotate("",xy=(x2,y2+.1),xytext=(x1+1.95,y1+.1),arrowprops=dict(arrowstyle="-|>",lw=.8,color="#838A92",shrinkA=5,shrinkB=5,connectionstyle="arc3,rad=0.03"))
    for enz,(x0,y0) in coords.items():
        ax.add_patch(FancyBboxPatch((x0,y0),2.15,.82,boxstyle="round,pad=.025,rounding_size=.08",fc="white",ec="#BAC2C9",lw=.75))
        ax.text(x0+.08,y0+.62,enz,fontsize=7.5,fontweight="bold",va="center")
        d=source[source.enzyme.eq(enz)].sort_values("stage_index")
        vals=np.vstack([d.RNA_z.to_numpy(),d.protein_z.to_numpy()])
        vals=np.nan_to_num(vals,nan=0)
        ax.imshow(vals,extent=(x0+.08,x0+2.05,y0+.08,y0+.43),aspect="auto",cmap=CMAP,norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2),interpolation="nearest",zorder=3)
        ax.text(x0+2.07,y0+.33,"R",fontsize=5.6,va="center",color=COL["red"])
        ax.text(x0+2.07,y0+.16,"P",fontsize=5.6,va="center",color=COL["purple"])
    ax.text(.55,.18,"Canonical biochemical topology; coloured strips are measured developmental profiles, not newly inferred causal regulation.",fontsize=7,color=COL["muted"])
    ax.text(10.05,1.15,"development",fontsize=6.8,color=COL["muted"])
    ax.imshow(np.linspace(-2,2,100).reshape(1,-1),extent=(11.25,14.55,1.12,1.34),aspect="auto",cmap=CMAP,norm=TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2))
    ax.text(11.25,.88,"low",ha="center",fontsize=6.5);ax.text(14.55,.88,"high",ha="center",fontsize=6.5)
    save3(fig,"Fig2g_canonical_lipid_pathway_multiomics_tracks")


def fig2h_candidates():
    cand=pd.read_csv(P["candidate"],sep="\t")
    rna=stage_mean_from_wide(pd.read_csv(P["rna_gene"],sep="\t"),["Geneid","Section","Enzyme","Preferred_name"])
    prot=stage_mean_from_wide(pd.read_csv(P["prot_gene"],sep="\t"),["Geneid","Section","Enzyme","Preferred_name"])
    common=set(rna.Geneid)&set(prot.Geneid)&set(cand.gene_id)
    chosen=cand[cand.gene_id.isin(common) & cand.protein_detected.eq(True) &
                pd.to_numeric(cand.current_protein_sample_completeness, errors="coerce").gt(0)].copy()
    chosen=chosen.sort_values(["figure2_evidence_count","current_protein_sample_completeness"],ascending=False).head(5)
    # If fewer than five pre-prioritised genes have both tracks, retain all available without fabrication.
    rr=rna.set_index("Geneid"); pp=prot.set_index("Geneid")
    rows=[]
    for _,r in chosen.iterrows():
        gid=r.gene_id
        rv=zscore(np.log2(rr.loc[gid,STAGES].astype(float).to_numpy()+1)); pv=zscore(np.log2(pp.loc[gid,STAGES].astype(float).to_numpy()+1))
        for i,s in enumerate(STAGES):
            rows.append({"gene_id":gid,"enzyme":rr.loc[gid,"Enzyme"],"stage":s,"stage_index":i+1,"RNA_z":rv[i],"protein_z":pv[i],
                         "RNA_module":r.rna_module,"RNA_kME":r.rna_module_kME,"protein_sample_completeness":r.current_protein_sample_completeness,
                         "selection_note":"pre-prioritised candidate with both gene-level RNA and protein tracks"})
    source=pd.DataFrame(rows); source.to_csv(OUT/"Fig2h_candidate_RNA_protein_trajectories_source.tsv",sep="\t",index=False)
    n=max(1,len(chosen)); fig,axes=plt.subplots(n,1,figsize=(8.2,1.1*n+1.45),sharex=True,sharey=True)
    if n==1:axes=[axes]
    for ax,(_,r) in zip(axes,chosen.iterrows()):
        d=source[source.gene_id.eq(r.gene_id)].sort_values("stage_index"); enz=d.enzyme.iloc[0]
        ax.plot(X,d.RNA_z,color=COL["red"],lw=1.8,marker="o",ms=2.6,label="RNA")
        ax.plot(X,d.protein_z,color=COL["purple"],lw=1.8,marker="s",ms=2.4,label="Protein")
        ax.axhline(0,color=COL["grid"],lw=.6)
        ax.text(.01,.81,f"{enz}  ·  {r.gene_id.replace('evm.TU.','')}",transform=ax.transAxes,fontweight="bold",fontsize=7.8)
        ax.text(.99,.81,f"RNA {r.rna_module}, kME={r.rna_module_kME:.2f}  |  protein completeness={r.current_protein_sample_completeness:.0%}",transform=ax.transAxes,ha="right",fontsize=6.4,color=COL["muted"])
        clean_ax(ax,None);ax.spines["bottom"].set_visible(False)
    axes[-1].spines["bottom"].set_visible(True);axes[-1].set_xticks(X,STAGES,rotation=45,ha="right");axes[-1].set_xlabel("Fruit developmental stage",labelpad=8)
    axes[len(axes)//2].set_ylabel("Within-track z-score")
    fig.suptitle("Representative candidates supported by measured RNA and protein trajectories",x=.09,ha="left",y=.985,fontsize=10,fontweight="bold")
    handles=[mpl.lines.Line2D([],[],color=COL["red"],marker="o",lw=1.8,label="RNA"),
             mpl.lines.Line2D([],[],color=COL["purple"],marker="s",lw=1.8,label="Protein")]
    fig.legend(handles=handles,frameon=False,ncol=2,loc="upper right",bbox_to_anchor=(.98,.987),fontsize=7)
    fig.text(.5,.018,"Candidate status is prioritisation, not functional validation. FAD2 mechanistic, ASE and haplotype analyses are reserved for Figure 3.",ha="center",fontsize=6.8,color=COL["muted"])
    panel_letter(axes[0],"h");fig.subplots_adjust(hspace=.08,bottom=.18,top=.91)
    save3(fig,"Fig2h_candidate_RNA_protein_trajectories")


def contact_sheet():
    stems=["Fig2a_compact_palm_phylogeny_CAFE5","Fig2b_CAFE5_GO_enrichment_conservative","Fig2c_core_lipid_gene_dosage_rigorous",
           "Fig2d_lipid_molecular_axes_13stage","Fig2e_lipid_RNA_developmental_programs","Fig2f_protein_lipid_execution_programme",
           "Fig2g_canonical_lipid_pathway_multiomics_tracks","Fig2h_candidate_RNA_protein_trajectories"]
    fig=plt.figure(figsize=(15,18));gs=fig.add_gridspec(4,2,height_ratios=[1,1.25,1,1.35],hspace=.08,wspace=.04)
    for i,stem in enumerate(stems):
        ax=fig.add_subplot(gs[i//2,i%2]);img=plt.imread(OUT/f"{stem}.png");ax.imshow(img);ax.axis("off")
    fig.suptitle("Figure 2 panel review sheet — evolution to developmental multi-omics",fontsize=15,fontweight="bold",y=.995)
    save3(fig,"Figure2_all_panels_review_sheet",dpi=250)


def write_manifest_and_notes(go_sets=None):
    rows=[]
    for label,path in P.items():
        p=Path(path);h=""
        if p.exists():
            sha=hashlib.sha256()
            with open(p,"rb") as fh:
                for block in iter(lambda:fh.read(1024*1024),b""):sha.update(block)
            h=sha.hexdigest()
        rows.append({"input":label,"path":str(p),"exists":p.exists(),"sha256":h})
    pd.DataFrame(rows).to_csv(OUT/"Figure2_input_manifest_sha256.tsv",sep="\t",index=False)
    notes="""# Figure 2 flat panel set: evidence and interpretation notes

## Main statement

Broad expansion of core lipid-gene dosage is not supported.  Instead, conserved
lipid genes show developmental RNA reprogramming accompanied by a lipid-enriched
protein programme and time-varying molecular lipid phenotypes.

## Panel-specific safeguards

- **2a:** CAFE5 numbers are branch changes, not branch-specific significant counts.
- **2b:** selected OGs are the intersection of oil-palm-ancestor branch changes and
  family-wide CAFE5 P<0.05.  The latter is not a branch-specific P value.
- **2c:** 107 unique lipid-associated OGs are deduplicated; five OGs map to two
  enzyme labels.  Audited outcome: 3 gained, 11 lost and 93 unchanged.
- **2d:** chemical evidence is shown explicitly. P02 and P06 are MS1-class proxies;
  named chemistry must not be interpreted as fully confirmed without MS2. The
  displayed trajectories pool two sampled genotypes (three replicates each;
  N=6 per stage).
- **2e:** only 13 fruit-development stages are clustered; postharvest stages were
  excluded from the high-oil developmental story.
- **2f:** the protein turquoise module is significantly enriched for lipid processes.
  No robust RNA-module–protein-module pair survived the prespecified permutation
  analysis, so no direct RNA→protein regulation is asserted.
- **2g:** arrows represent canonical lipid biochemistry, while coloured strips are
  measured RNA/protein trajectories.  This is a data overlay, not a new mechanism.
- **2h:** candidates are prioritised, not functionally validated. FAD2 ASE/haplotype
  interpretation belongs in Figure 3.

## File structure

Every PDF, SVG, PNG, source table, manifest and script is directly in this folder.
No nested output folder is created.
"""
    (OUT/"README_Figure2_panels_and_evidence.md").write_text(notes,encoding="utf-8")


def main():
    setup_style()
    write_manifest_and_notes()
    fig2a_tree()
    go_df,go_sets,bg=go_enrichment()
    subprocess.run(["Rscript", str(OUT / "render_Fig2b_GO_enrichment.R")], check=True)
    fig2c_copy_number();fig2d_axes();fig2e_rna_programs();fig2f_protein();fig2g_pathway();fig2h_candidates()
    contact_sheet();write_manifest_and_notes(go_sets)
    print(f"Rendered Figure 2 panels into: {OUT}")


if __name__=="__main__":
    main()
