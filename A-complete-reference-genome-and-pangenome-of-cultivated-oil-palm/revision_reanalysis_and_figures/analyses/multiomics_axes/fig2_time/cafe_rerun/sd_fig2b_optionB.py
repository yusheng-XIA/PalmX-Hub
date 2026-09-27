"""Source Data Fig2b_gene_families from a CAFE5 gamma run, in the layout of the current sheet (node ids of the
authors' run, i.e. the ids printed in the current sheet) plus clade and MCMCTree node columns."""
import re, sys
run = sys.argv[1]
ORDER = [(0, "Pisifera", {"Pisifera"}), (1, "Dura", {"Dura"}), (2, "", {"Dura", "Pisifera"}), (3, "American_hap1", {"American_hap1"}),
         (4, "", {"Dura", "Pisifera", "American_hap1"}), (5, "Cocos_nucifera", {"Cocos_nucifera"}),
         (6, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera"}), (7, "Areca_catechu", {"Areca_catechu"}),
         (8, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu"}), (9, "Phoenix_dactylifera", {"Phoenix_dactylifera"}),
         (10, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera"}),
         (11, "Nypa_fruticans", {"Nypa_fruticans"}), (12, "Daemonorops", {"Daemonorops"}), (13, "Calamus", {"Calamus"}),
         (14, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans"}),
         (15, "", {"Calamus", "Daemonorops"}), (16, "Musa_balbisiana", {"Musa_balbisiana"}), (17, "Musa_acuminata", {"Musa_acuminata"}),
         (18, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans", "Calamus", "Daemonorops"}),
         (19, "", {"Musa_acuminata", "Musa_balbisiana"}), (20, "", {"Dura", "Pisifera", "American_hap1", "Cocos_nucifera", "Areca_catechu", "Phoenix_dactylifera", "Nypa_fruticans", "Calamus", "Daemonorops", "Musa_acuminata", "Musa_balbisiana"}),
         (21, "Oryza_sativa", {"Oryza_sativa"}), (22, "", None)]
TN = {2: "t_n43", 4: "t_n42", 6: "t_n41", 8: "t_n40", 10: "t_n39", 14: "t_n38", 15: "t_n37", 18: "t_n36", 19: "t_n35", 20: "t_n34", 22: "t_n33"}
CL = {2: "Dura + pisifera", 4: "Elaeis (FL-Hap1 + dura + pisifera)", 6: "Cocos + Elaeis", 8: "Areca + (Cocos, Elaeis)",
      10: "Phoenix + (Areca, Cocos, Elaeis)", 14: "Nypa + other palms", 15: "Calamus + Daemonorops", 18: "Palm crown", 19: "Musa",
      20: "Musa + palms", 22: "Root (Oryza + Musa + palms)"}
line = next(l for l in open(run + "/Gamma_asr.tre") if "TREE" in l and "=" in l)
t = line.split("=", 1)[1]
out, cur = {}, [[]]
for m in re.finditer(r"\(|\)<(\d+)>|([A-Za-z][A-Za-z_0-9]*)<(\d+)>", t):
    s = m.group(0)
    if s == "(": cur.append([])
    elif s.startswith(")"):
        k = cur.pop(); cur[-1].extend(k); out[frozenset(k)] = m.group(1)
    else: cur[-1].append(m.group(2)); out[frozenset([m.group(2)])] = m.group(3)
root = max(out, key=len); out_root = out[root]
hdr, *rows = [l.rstrip("\n").split("\t") for l in open(run + "/Gamma_change.tab")]
print("node_id\tcolumn\tclade\tMCMCTree node (final run)\texpanded_families\tcontracted_families\tmaintained_families")
for aid, tip, s in ORDER:
    rid = out_root if s is None else out[frozenset(s)]
    j = next(i for i, h in enumerate(hdr) if h.endswith(f"<{rid}>"))
    v = [int(r[j]) for r in rows]
    print(f"{aid}\t{tip}<{aid}>\t{CL.get(aid, tip.replace('_', ' '))}\t{TN.get(aid, '')}\t{sum(x > 0 for x in v)}\t{sum(x < 0 for x in v)}\t{sum(x == 0 for x in v)}")
