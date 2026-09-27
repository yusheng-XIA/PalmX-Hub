"""SF5 (sfplan): d, label the OLE16b tip whose protein is identical in FL-Hap2 and TN-Hap2; i, number of OLE16a-positive FL nuclei above each bar."""
import fitz, sys
src, dst = sys.argv[1], sys.argv[2]
FR = "/System/Library/Fonts/Supplemental/Arial.ttf"; FB = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
d = fitz.open(src); p = d[0]
# --- d: OLE16b FL-Hap2 (chr04B) tip = identical protein in TN-Hap2 (unified UFTN008332: evm.TU.chr04B.697; bk_hap2_chr3.693)
spans = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
t = [s for s in spans if s["text"] == "OLE16b FL-Hap2 (chr04B)"]; assert len(t) == 1
s = t[0]; x, y = s["bbox"][2], s["origin"][1]
p.insert_font(fontname="AB", fontfile=FB); p.insert_font(fontname="AR", fontfile=FR)
p.insert_text((x, y), " = TN-Hap2", fontname="AB", fontsize=5, color=(0x55/255,) * 3)
# --- i: positive-nucleus counts (Source Data SF5i: % x nuclei), FL 185 d
cnt = {"C18": 24, "C12": 5, "C7": 79, "C15": 8, "C9": 56, "C4": 3, "C0": 5, "C14": 2, "C17": 1, "C5": 1, "C13": 1,
       "C6": 1, "C1": 2, "C2": 0, "C3": 0, "C8": 0, "C16": 0}
assert sum(cnt.values()) == 188
labs = {s["text"]: (s["bbox"][0] + s["bbox"][2]) / 2 for s in spans if s["text"] in cnt and s["bbox"][1] > 620 and s["bbox"][0] > 380}
assert len(labs) == 17
teal = (0.137, 0.537, 0.49)
bars = [dr["rect"] for dr in p.get_drawings() if dr.get("fill") and all(abs(a - b) < 0.01 for a, b in zip(dr["fill"], teal))
        and 520 < dr["rect"].y0 < 640 and 395 < dr["rect"].x0 < 500 and dr["rect"].width < 3]
assert len(bars) == 17, len(bars)
for cl, xc in labs.items():
    r = min(bars, key=lambda r: abs((r.x0 + r.x1) / 2 - xc))
    top = r.y0 if r.height > 0.05 else 628.0
    p.insert_text((r.x1 - 0.2, top - 1.0), str(cnt[cl]), fontname="AR", fontsize=5, rotate=90, color=(0.2, 0.2, 0.2))
d.save(dst, garbage=3, deflate=True)
