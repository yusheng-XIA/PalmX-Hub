"""Fig. 3/5 minor in-figure edits (vector; text-only redaction, no fill, graphics untouched).
Fig. 3: 'Postharvest' -> 'Post-harvest' (3 spans: 3c header, 3h header, 3i legend), same baseline/size/colour/weight,
        re-centred (3c, 3h) or left-aligned (3i).  'Same CDS' (3a) kept: verified as identical CDS (see README).
Fig. 5: remove the 47 red dSV/dSNP labels in 5i (red boxes kept; chr01B totals remain in 5h, per-window counts in
        Source Data Fig5i_selected_chr01B).  Outputs fix/fig35_minor/Figure{3,5}.pdf/.tif + before/after crops."""
import fitz, io, re
from PIL import Image
S = "deliver/Main_Figures_revised/"; O = "fix/fig35_minor/"
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"; ARB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
norm = lambda t: t.replace("\xa0", " ")

def spans(p):
    return [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]

def rgb(c): return ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255

def redact(p, ss):
    for s in ss:
        x0, y0, x1, y1 = s["bbox"]; h = y1 - y0
        p.add_redact_annot(fitz.Rect(x0 + 0.1, y0 + 0.3 * h, x1 - 0.1, y1 - 0.3 * h), fill=False)
    p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                       text=fitz.PDF_REDACT_TEXT_REMOVE)

def export(doc, name, crops):
    doc.save(O + name + ".pdf", garbage=3, deflate=True)
    pg = fitz.open(O + name + ".pdf")[0]
    pm = pg.get_pixmap(matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False)
    Image.open(io.BytesIO(pm.tobytes("png"))).convert("RGB").save(O + name + ".tif", compression="tiff_lzw", dpi=(600, 600))
    before = fitz.open(S + name + ".pdf")[0]
    for tag, r in crops.items():
        a = before.get_pixmap(clip=fitz.Rect(r), dpi=400); b = pg.get_pixmap(clip=fitz.Rect(r), dpi=400)
        A = Image.open(io.BytesIO(a.tobytes("png"))); B = Image.open(io.BytesIO(b.tobytes("png")))
        C = Image.new("RGB", (A.width, A.height * 2 + 20), "white"); C.paste(A, (0, 0)); C.paste(B, (0, A.height + 20))
        C.save(O + f"{name}_{tag}_before_after.png")
    return pg, pm

# ---------- Figure 3 ----------
d = fitz.open(S + "Figure3.pdf"); p = d[0]
tgt = [s for s in spans(p) if "Postharvest" in norm(s["text"])]
assert len(tgt) == 3, [s["text"] for s in tgt]
plan = []
for s in tgt:
    new = norm(s["text"]).replace("Postharvest", "Post-harvest").strip()
    bold = "Bold" in s["font"]
    align = "left" if s["bbox"][1] > 390 else "center"   # 3i legend entry is left-aligned
    plan.append((s, new, bold, align))
redact(p, tgt)
for k, (s, new, bold, align) in enumerate(plan):
    ff = ARB if bold else AR
    w = fitz.Font(fontfile=ff).text_length(new, fontsize=s["size"])
    x0, y0, x1, y1 = s["bbox"]
    x = x0 if align == "left" else (x0 + x1) / 2 - w / 2
    p.insert_text((x, s["origin"][1]), new, fontname=f"Fm35{k}", fontfile=ff, fontsize=s["size"], color=rgb(s["color"]))
pg3, pm3 = export(d, "Figure3", {"3c": (300, 18, 424, 40), "3h": (330, 240, 424, 258), "3i": (20, 385, 150, 410)})
t_old = norm(fitz.open(S + "Figure3.pdf")[0].get_text()).replace("\xad", "-"); t_new = norm(pg3.get_text()).replace("\xad", "-")
assert sorted(re.sub("Postharvest", "Post-harvest", t_old).split()) == sorted(t_new.split()), "Fig3 other text changed"
assert "Postharvest" not in t_new and t_new.count("Post-harvest") == 3
print("Figure3 ok", pm3.width, pm3.height)

# ---------- Figure 5 ----------
d = fitz.open(S + "Figure5.pdf"); p = d[0]
lab = [s for s in spans(p) if s["color"] == 0xb83f3e and re.fullmatch(r"\d+/\d+", s["text"].strip()) and s["bbox"][1] > 455 and s["bbox"][0] > 280]
assert len(lab) == 47, len(lab)
dsv = sum(int(s["text"].split("/")[0]) for s in lab); dsnp = sum(int(s["text"].split("/")[1]) for s in lab)
print("removed labels sum dSV/dSNP =", dsv, dsnp)
redact(p, lab)
pg5, pm5 = export(d, "Figure5", {"5i": (240, 455, 518, 695)})
o = [norm(s["text"]) for s in spans(fitz.open(S + "Figure5.pdf")[0])]
n = [norm(s["text"]) for s in spans(pg5)]
removed = [s["text"] for s in lab]
for t in removed: o.remove(t)
assert sorted(o) == sorted(n), "Fig5 other text changed"
print("Figure5 ok", pm5.width, pm5.height)
