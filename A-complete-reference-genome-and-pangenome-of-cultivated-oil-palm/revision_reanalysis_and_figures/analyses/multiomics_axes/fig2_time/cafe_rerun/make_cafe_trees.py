"""Build 12-taxon integer-Ma ultrametric CAFE5 trees (same topology and tip set as the authors' layer2_palm
cafe_tree.nwk) from the old (2026-02-23) and new (2026-06-17) MCMCTree FigTree.tre posterior-mean trees.
Node ages are x100 (MCMCTree unit = 100 Ma) and rounded to integers; branch = parent age - child age."""
import re, sys
TAXA = ["Oryza_sativa","Musa_acuminata","Musa_balbisiana","Calamus","Daemonorops","Nypa_fruticans",
        "Phoenix_dactylifera","Areca_catechu","Cocos_nucifera","American_hap1","Dura","Pisifera"]
TOPO = "(Oryza_sativa,((Musa_acuminata,Musa_balbisiana),((Calamus,Daemonorops),(Nypa_fruticans,(Phoenix_dactylifera,(Areca_catechu,(Cocos_nucifera,(American_hap1,(Dura,Pisifera)))))))));"

def parse(s):
    s = re.sub(r"\[.*?\]", "", s).replace(" ", "")
    pos = 0
    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1; kids = [node()]
            while s[pos] == ",": pos += 1; kids.append(node())
            pos += 1
            name = ""
        else:
            m = re.match(r"[A-Za-z_0-9]+", s[pos:]); name = m.group(0); pos += len(name); kids = []
        m = re.match(r"(?:[A-Za-z_0-9]*)?:([0-9.eE+-]+)", s[pos:])
        bl = 0.0
        if m: bl = float(m.group(1)); pos += len(m.group(0))
        return {"name": name, "kids": kids, "bl": bl}
    return node()

def tips(n): return [n["name"]] if not n["kids"] else sum((tips(k) for k in n["kids"]), [])
def height(n): return 0.0 if not n["kids"] else height(n["kids"][0]) + n["kids"][0]["bl"]

def mrca_age(root, group):
    best = None
    def walk(n):
        nonlocal best
        t = set(tips(n))
        if set(group) <= t:
            best = n
            for k in n["kids"]: walk(k)
    walk(root)
    return height(best)

def build(tree_text):
    line = next(l for l in tree_text.splitlines() if "UTREE" in l)
    root = parse(line.split("=", 1)[1].strip().rstrip(";"))
    topo = parse(TOPO.rstrip(";"))
    ages = {}
    def rec(n):
        if not n["kids"]: return n["name"], 0
        parts = [rec(k) for k in n["kids"]]
        a = round(100 * mrca_age(root, tips(n)))
        return "(" + ",".join(f"{p}:{a - pa}" for p, pa in parts) + ")", a
    s, a = rec(topo)
    # node table
    tab = []
    def rec2(n):
        if n["kids"]:
            tab.append(("+".join(tips(n)[:2]) + ("..." if len(tips(n)) > 2 else ""), 100 * mrca_age(root, tips(n))))
            for k in n["kids"]: rec2(k)
    rec2(topo)
    return s + ";", tab

for tag in ("old", "new"):
    nwk, tab = build(open(f"{tag}_FigTree.tre").read())
    open(f"cafe_tree_{tag}MCMC.nwk", "w").write(nwk + "\n")
    print(tag, nwk)
    for c, a in tab: print(f"  {c}\t{a:.2f}")
