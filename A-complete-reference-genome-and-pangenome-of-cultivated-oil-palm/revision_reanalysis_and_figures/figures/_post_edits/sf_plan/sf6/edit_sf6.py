"""SF6 (sfplan S26): conditioning-QC squares in b were hidden beneath the analysis-QC triangles; overlay dark open squares
(same positions, 0.6 pt larger) so the three conditioning injections are visible; legend symbol styled the same."""
import sys, fitz
src, dst = sys.argv[1], sys.argv[2]
d = fitz.open(src); p = d[0]
sq = [dr["rect"] for dr in p.get_drawings() if dr.get("fill") and tuple(round(c, 2) for c in dr["fill"]) == (0.73, 0.75, 0.78)
      and dr["rect"].width < 6]
assert len(sq) == 7, len(sq)
sh = p.new_shape()
for r in sq:
    rr = fitz.Rect(r.x0 - 0.6, r.y0 - 0.6, r.x1 + 0.6, r.y1 + 0.6)
    sh.draw_rect(rr)
sh.finish(color=(0.25, 0.27, 0.30), fill=None, width=0.45)
sh.commit(overlay=True)
d.save(dst, garbage=3, deflate=True); print("overlaid", len(sq))
