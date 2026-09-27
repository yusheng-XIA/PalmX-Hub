#!/usr/bin/env python3
"""Before/after comparison, vector-text diff and deterministic checks for the beautify pass.

python3 figtools.py compare OLD.pdf NEW.pdf KEY     -> compare/KEY_before_after.png, compare/KEY_textdiff.tsv
python3 figtools.py check NEW.pdf NEW.png KEY TYPE  -> checks/KEY_check.txt (+ critic json)
"""
import collections, json, re, subprocess, sys
from pathlib import Path

import fitz
from PIL import Image, ImageDraw, ImageFont

B = Path(__file__).resolve().parents[1]
SKILL = B.parents[1] / "skill_make_figures/make-figures/scripts/critic_figure.py"
MM = 25.4 / 72
SUP = set("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺")


def spans(pdf):
    out = []
    pg = fitz.open(pdf)[0]
    for b in pg.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip():
                    out.append(s)
    return out, pg


def words(pdf):
    pg = fitz.open(pdf)[0]
    return [w[4] for w in pg.get_text("words")]


def norm_tokens(pdf):
    """Word tokens with glyph-level normalisation that does not change meaning."""
    toks = []
    for w in words(pdf):
        w = w.replace("−", "-").replace(" ", " ")
        toks += [t for t in re.split(r"\s+", w) if t]
    return collections.Counter(toks)


def render(pdf, dpi):
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False)
    return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)


def compare(old, new, key, dpi=130):
    a, b = render(old, dpi), render(new, dpi)
    h = max(a.height, b.height)
    gap, top = 30, 46
    im = Image.new("RGB", (a.width + b.width + gap, h + top), "white")
    im.paste(a, (0, top)); im.paste(b, (a.width + gap, top))
    d = ImageDraw.Draw(im)
    try:
        f = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 26)
    except OSError:
        f = None
    d.text((10, 8), f"{key}  BEFORE (current submission)", fill="black", font=f)
    d.text((a.width + gap + 10, 8), f"{key}  AFTER (beautify candidate)", fill="black", font=f)
    d.line([(a.width + gap // 2, 0), (a.width + gap // 2, h + top)], fill=(150, 150, 150), width=2)
    (B / "compare").mkdir(exist_ok=True)
    im.save(B / "compare" / f"{key}_before_after.png", optimize=True)
    # text diff
    o, n = norm_tokens(old), norm_tokens(new)
    rem, add = o - n, n - o
    rows = [("removed(before_only)", t, c) for t, c in sorted(rem.items())] + \
           [("added(after_only)", t, c) for t, c in sorted(add.items())]
    num = lambda t: bool(re.search(r"\d", t))
    with open(B / "compare" / f"{key}_textdiff.tsv", "w") as fh:
        fh.write(f"# tokens before={sum(o.values())} after={sum(n.values())}; "
                 f"numeric tokens changed={sum(c for s, t, c in rows if num(t))}\n")
        fh.write("side\ttoken\tcount\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{key}: tokens {sum(o.values())}->{sum(n.values())}; removed {sum(rem.values())}, added {sum(add.values())}; "
          f"numeric diffs: {[r for r in rows if num(r[1])][:12]}")


def check(pdf, png, key, typ, width_in=7.087):
    sp, pg = spans(pdf)
    W, H = pg.rect.width * MM, pg.rect.height * MM
    small = [(round(s["size"], 2), s["text"].strip()) for s in sp
             if s["size"] < 4.95 and not (set(s["text"].strip()) <= SUP | set("0123456789−-+ ."))]
    small_sup = [(round(s["size"], 2), s["text"].strip()) for s in sp if s["size"] < 4.95]
    letters = [(s["text"].strip(), round(s["size"], 2), s["font"]) for s in sp
               if re.fullmatch(r"[a-z]", s["text"].strip()) and "Bold" in s["font"] and s["size"] >= 7.5]
    fonts = collections.Counter(s["font"] for s in sp)
    thin = collections.Counter()
    for d in pg.get_drawings():
        w = d.get("width")
        if d.get("color") is not None and w is not None and 0 < w < 0.25:
            thin[round(w, 3)] += 1
    rep = [f"# {key}: {W:.1f} x {H:.1f} mm; spans {len(sp)}",
           f"fonts: {dict(fonts)}",
           f"panel letters: {letters}",
           f"body text < 5 pt (excluding pure super/subscript digits): {len(small)} {small[:25]}",
           f"all spans < 5 pt incl. super/subscripts: {len(small_sup)} sizes {sorted(set(s for s, _ in small_sup))}",
           f"stroked paths with width < 0.25 pt: {dict(thin)}"]
    (B / "checks").mkdir(exist_ok=True)
    js = B / "checks" / f"{key}_critic.json"
    r = subprocess.run([sys.executable, str(SKILL), str(png), "--type", typ, "--spec-width-in", str(width_in),
                        "--exploratory", "--out", str(js)] + ([] if typ == "heatmap" else ["--art-class", "line_art"]), capture_output=True, text=True)
    try:
        j = json.loads(js.read_text())
        rep.append(f"critic_figure.py --exploratory: {j.get('status')} | {j.get('summary')}")
        for f in j.get("flags", []):
            rep.append(f"  flag: {f}")
    except Exception as e:  # noqa
        rep.append("critic: " + (r.stdout[-800:] + r.stderr[-800:]))
    (B / "checks" / f"{key}_check.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    if sys.argv[1] == "compare":
        compare(*sys.argv[2:5])
    else:
        check(*sys.argv[2:6])


def docx_png(png, out):
    """2400-px-wide RGB copy (same convention as docx_img/)."""
    im = Image.open(png)
    if im.mode in ("RGBA", "LA"):
        bg = Image.new("RGB", im.size, "white"); bg.paste(im, mask=im.split()[-1]); im = bg
    else:
        im = im.convert("RGB")
    h = round(im.height * 2400 / im.width)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    im.resize((2400, h), Image.LANCZOS).save(out, optimize=True)
