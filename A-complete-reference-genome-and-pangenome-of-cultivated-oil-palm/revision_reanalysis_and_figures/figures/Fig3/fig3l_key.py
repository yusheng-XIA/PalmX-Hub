"""Item 11: Fig. 3l key header 'Favourable action (293 targets)' also sits over the magenta 'Haplotype screen' diamond,
which marks the 43 actual-haplotype screens (legend l; Source Data Fig5h_Fig3l_targets: 145 introduce/tune + 128 retain FL +
20 timing screen = 293 triangles; 43 'Validate only' actual-haplotype targets = diamonds).  Header -> 'Target action (336 targets)',
re-centred on the original header.  Only this span changes.  Output fix/orphan_fix/Figure3.pdf/.tif (600 dpi, LZW, RGB)."""
import fitz, io
from PIL import Image
SRC = "deliver/Main_Figures_revised/Figure3.pdf"; OUT = "fix/orphan_fix/Figure3"
OLD, NEW = "Favourable action (293 targets)", "Target action (336 targets)"
AR = "/System/Library/Fonts/Supplemental/Arial.ttf"
d = fitz.open(SRC); p = d[0]
norm = lambda t: t.replace("\xa0", " ")
hit = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if norm(s["text"]) == OLD]
assert len(hit) == 1, len(hit)
s = hit[0]; r = fitz.Rect(s["bbox"])
p.add_redact_annot(r, fill=(1, 1, 1)); p.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE)
font = fitz.Font(fontfile=AR)
w_new = font.text_length(NEW, fontsize=s["size"])
cx = (r.x0 + r.x1) / 2
c = s["color"]; rgb = ((c >> 16) & 255) / 255, ((c >> 8) & 255) / 255, (c & 255) / 255
p.insert_text((cx - w_new / 2, s["origin"][1]), NEW, fontname="Fol1", fontfile=AR, fontsize=s["size"], color=rgb)
d.save(OUT + ".pdf", garbage=3, deflate=True)
pg = fitz.open(OUT + ".pdf")[0]
assert pg.search_for("Target action (336 targets)") and not pg.search_for("Favourable action")
old = fitz.open(SRC)[0].get_text().replace("\xa0", " ").replace(OLD, "")
new = pg.get_text().replace("\xa0", " ").replace(NEW, "")
assert old.split() == new.split(), "other text changed"
pm = pg.get_pixmap(matrix=fitz.Matrix(600 / 72, 600 / 72), alpha=False)
Image.open(io.BytesIO(pm.tobytes("png"))).convert("RGB").save(OUT + ".tif", compression="tiff_lzw", dpi=(600, 600))
pg.get_pixmap(clip=fitz.Rect(250, 570, 424, 625), dpi=600).save("fix/orphan_fix/Figure3l_key_after.png")
print("ok", pm.width, pm.height)
