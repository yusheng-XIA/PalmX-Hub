"""Compare CAFE5 gamma results: author run (layer2_palm/gamma_results) vs re-runs on the original tree, the old-run
MCMCTree tree and the final-run MCMCTree tree. Node IDs are mapped by clade (CAFE 5.1 numbers nodes differently)."""
import re, sys, csv
A = "${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder/phylo_divtime_cafe/03_cafe5/layer2_palm/gamma_results"
R = "${CLUSTER_WORK}/fig2_time/runs"
RUNS = [("author_orig_tree", A), ("rerun_orig_tree", R + "/gamma_orig"), ("rerun_oldMCMC_tree", R + "/gamma_oldMCMC"),
        ("rerun_newMCMC_tree", R + "/gamma_newMCMC")]
def clades(asr):
    line = next(l for l in open(asr) if "TREE" in l and "=" in l)
    t = line.split("=", 1)[1].strip()
    stack, out = [], {}
    tok = re.finditer(r"\(|\)<(\d+)>|([A-Za-z][A-Za-z_0-9]*)<(\d+)>|,|;", t)
    cur = [[]]
    for m in tok:
        s = m.group(0)
        if s == "(": cur.append([])
        elif s.startswith(")"):
            kids = cur.pop(); cur[-1].extend(kids); out[m.group(1)] = frozenset(kids)
        elif m.group(2): cur[-1].append(m.group(2)); out[m.group(3)] = frozenset([m.group(2)])
    return {v: k for k, v in out.items()}
NAMES = [("DP", {"Dura", "Pisifera"}), ("Elaeis", {"Dura", "Pisifera", "American_hap1"}),
         ("Cocos+Elaeis", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera"}),
         ("Areca+", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu"}),
         ("Phoenix+", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera"}),
         ("Nypa+", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans"}),
         ("Calamus+Daemonorops", {"Calamus", "Daemonorops"}),
         ("Palm crown", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans", "Calamus", "Daemonorops"}),
         ("Musa", {"Musa_acuminata", "Musa_balbisiana"}),
         ("Musa+palms", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans", "Calamus", "Daemonorops", "Musa_acuminata", "Musa_balbisiana"})]
TIPS = ["Pisifera", "Dura", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans",
        "Daemonorops", "Calamus", "Musa_balbisiana", "Musa_acuminata", "Oryza_sativa"]
rows = {}; meta = {}
for name, d in RUNS:
    try:
        cmap = clades(d + "/Gamma_asr.tre")
    except (FileNotFoundError, StopIteration):
        continue
    cr = {}
    for l in open(d + "/Gamma_clade_results.txt"):
        if l.startswith("#"): continue
        k, inc, dec = l.rstrip("\n").split("\t"); cr[re.search(r"<(\d+)>", k).group(1)] = (int(inc), int(dec))
    for lab, s in NAMES + [(t, {t}) for t in TIPS]:
        rows.setdefault(lab, {})[name] = cr[cmap[frozenset(s)]]
    res = open(d + "/Gamma_results.txt").read()
    lam = float(re.search(r"Lambda: ([\d.eE-]+)", res).group(1)); lnl = re.search(r"-lnL\): ([\d.]+)", res).group(1)
    alpha = re.search(r"Alpha: ([\d.eE-]+)", res)
    nsig = sum(1 for l in open(d + "/Gamma_family_results.txt") if not l.startswith("#") and l.split("\t")[1] not in ("", "pvalue")
               and float(l.split("\t")[1]) < 0.05)
    meta[name] = (lam, alpha.group(1) if alpha else "", lnl, nsig)
w = csv.writer(sys.stdout, delimiter="\t")
names = [n for n, _ in RUNS if n in meta]
w.writerow(["node"] + [f"{n}_{x}" for n in names for x in ("expanded", "contracted")])
for lab in [n for n, _ in NAMES] + TIPS:
    w.writerow([lab] + [v for n in names for v in rows[lab][n]])
for i, k in enumerate(["lambda", "alpha", "-lnL", "families_P<0.05"]):
    w.writerow([k] + [x for n in names for x in (meta[n][i], "")])
