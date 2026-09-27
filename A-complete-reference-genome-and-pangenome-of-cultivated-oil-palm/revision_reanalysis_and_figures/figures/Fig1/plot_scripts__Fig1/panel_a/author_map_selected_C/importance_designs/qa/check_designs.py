"""600 dpi geometry/data checks; run from fig1a_draft/."""
from pathlib import Path
import json
import re
import numpy as np
from matplotlib.text import Text, Annotation
from matplotlib.patches import Rectangle, ConnectionPatch
from shapely.geometry import box
from PIL import Image

p = Path("plot_oilpalm_importance.py")
ns = {}
exec(compile(re.sub(r"fig\.savefig\([^\n]*\)", "pass", p.read_text()), str(p), "exec"), ns)
reports = {}
for key, fig in ns["figures"].items():
    fig.set_dpi(600); fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    texts = list(fig.texts)
    for ax in fig.axes:
        texts += list(ax.texts) + ax.get_xticklabels() + ax.get_yticklabels() + [ax.xaxis.label, ax.yaxis.label, ax.title]
        if ax.get_legend(): texts += list(ax.get_legend().get_texts())
    for legend in fig.legends: texts += list(legend.get_texts())
    texts = [t for t in texts if t.get_visible() and t.get_text().strip()]
    boxes = [Text.get_window_extent(t, renderer=renderer) for t in texts]
    flags = []
    for i, (text, bb) in enumerate(zip(texts, boxes)):
        if not (0 <= bb.x0 <= bb.x1 <= fig.bbox.width and 0 <= bb.y0 <= bb.y1 <= fig.bbox.height):
            flags.append("outside figure: " + text.get_text())
        if text.get_fontsize() < 5.5: flags.append("font below 5.5 pt")
        for j in range(i+1, len(texts)):
            if bb.overlaps(boxes[j]): flags.append(f"text overlap: {text.get_text()} / {texts[j].get_text()}")
        for other in texts:
            if other is not text and other.get_bbox_patch() is not None:
                if bb.overlaps(other.get_bbox_patch().get_window_extent(renderer)):
                    flags.append("text on badge: " + text.get_text())
        if text.axes:
            for patch in text.axes.patches:
                if isinstance(patch, Rectangle) and patch.get_zorder()>=8 and bb.overlaps(patch.get_window_extent(renderer)):
                    flags.append("text on scale bar: " + text.get_text())
    maps = [a for a in fig.axes if tuple(round(v, 3) for v in a.get_xlim()) in
            [(-118.,168.),(-104.,-44.),(-16.,32.),(92.,156.)]]
    for ax in maps:
        for text in ax.texts:
            if not text.get_text(): continue
            bb = Text.get_window_extent(text, renderer=renderer)
            corners = ax.transData.inverted().transform([[bb.x0,bb.y0],[bb.x1,bb.y1]])
            if any(ns["map_ns"]["land_all"].intersects(box(*corners.flatten()))):
                flags.append("map label on land: " + text.get_text())
            for other in ax.texts:
                if other is text or not isinstance(other, Annotation) or other.arrow_patch is None: continue
                patch = other.arrow_patch
                path = patch.get_path().transformed(patch.get_transform())
                if bb.overlaps(patch.get_window_extent(renderer)) and path.intersects_bbox(bb, filled=False):
                    flags.append(f"leader on text: {other.get_text()} / {text.get_text()}")
    for artist in fig.artists:
        if isinstance(artist, ConnectionPatch):
            path = artist.get_path().transformed(artist.get_transform())
            for text, bb in zip(texts, boxes):
                if path.intersects_bbox(bb, filled=False): flags.append("zoom line on text: " + text.get_text())
    im = Image.open(Path("importance_designs") / f"Fig_oilpalm_importance_{key}.png")
    reports[key] = {"size_mm": (fig.get_size_inches()*25.4).tolist(),
                    "pixels": im.size, "dpi": im.info.get("dpi"), "text_count":len(texts),
                    "minimum_font_pt":min(t.get_fontsize() for t in texts), "findings":flags}
    assert abs(fig.get_figwidth()*25.4-183)<1e-6 and fig.get_figheight()*25.4<=80
    fig.clear()
summary = json.loads(Path("importance_designs/data_summary.json").read_text())
assert len(summary["area_categories"])==9 and len(summary["oil_categories"])==11
assert np.isclose(summary["oil_palm_area_share_pct"], 100*summary["oil_palm_area_mha"]/summary["area_total_mha"])
assert np.isclose(summary["oil_palm_oil_share_pct"], 100*summary["oil_palm_oil_mt"]/summary["oil_total_mt"])
Path("importance_designs/qa/geometry_600dpi.json").write_text(json.dumps(reports,ensure_ascii=False,indent=2))
print(json.dumps(reports,ensure_ascii=False,indent=2))
if any(r["findings"] for r in reports.values()): raise SystemExit(1)
