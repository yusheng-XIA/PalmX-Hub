#!/usr/bin/env python3
"""Ultrametric CAFE5 input tree for the 12 monocot genomes, pruned from the MCMCTree posterior-mean tree.

Node ages (MCMCTree unit = 100 Ma) are converted to Ma and rounded to whole million years;
branch length = parent age - child age.
usage: 03_make_cafe_tree.py FigTree.tre > cafe_tree.nwk
"""
import re
import sys

TOPO = ("(Oryza_sativa,((Musa_acuminata,Musa_balbisiana),((Calamus,Daemonorops),(Nypa_fruticans,"
        "(Phoenix_dactylifera,(Areca_catechu,(Cocos_nucifera,(American_hap1,(Dura,Pisifera)))))))));")
# American_hap1 = FL-Hap1, Dura = TK, Pisifera = NS


def parse(s):
    s = re.sub(r"\[.*?\]", "", s).replace(" ", "")
    pos = 0

    def node():
        nonlocal pos
        kids, name = [], ""
        if s[pos] == "(":
            pos += 1; kids.append(node())
            while s[pos] == ",":
                pos += 1; kids.append(node())
            pos += 1
        else:
            name = re.match(r"[A-Za-z_0-9]+", s[pos:]).group(0); pos += len(name)
        m = re.match(r"(?:[A-Za-z_0-9]*)?:([0-9.eE+-]+)", s[pos:])
        bl = 0.0
        if m:
            bl = float(m.group(1)); pos += len(m.group(0))
        return {"name": name, "kids": kids, "bl": bl}
    return node()


def tips(n):
    return [n["name"]] if not n["kids"] else sum((tips(k) for k in n["kids"]), [])


def height(n):
    return 0.0 if not n["kids"] else height(n["kids"][0]) + n["kids"][0]["bl"]


def mrca_age(root, group):
    group = set(group)

    def walk(n):
        for k in n["kids"]:
            if group <= set(tips(k)):
                return walk(k)
        return n
    return height(walk(root))


def build(tree_text):
    full = parse(tree_text.strip().split("=", 1)[-1].rstrip(";"))
    topo = parse(TOPO.rstrip(";"))

    def age(n):
        return 0 if not n["kids"] else int(round(100 * mrca_age(full, tips(n))))

    def emit(n, parent_age=None):
        a = age(n)
        inner = n["name"] if not n["kids"] else "(" + ",".join(emit(k, a) for k in n["kids"]) + ")"
        return inner if parent_age is None else f"{inner}:{parent_age - a}"
    return emit(topo) + ";"


if __name__ == "__main__":
    text = open(sys.argv[1]).read()
    tree_line = next(l for l in text.splitlines() if "(" in l)
    print(build(tree_line))
