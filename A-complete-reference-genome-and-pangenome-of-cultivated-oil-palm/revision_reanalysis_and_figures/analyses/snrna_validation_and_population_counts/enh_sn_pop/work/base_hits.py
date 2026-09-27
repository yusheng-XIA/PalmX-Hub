import sys, re, pickle
from pathlib import Path
H = Path(__file__).resolve().parent; S = H.parents[2]
sys.path[:0] = [str(H / "pb")]
import build_main_text as BM
from docx_accept import accept_all
d = BM.Doc(S / "v1639/HB-01_MS_FinalFigures_HighRes_Legends_Tracked.docx"); accept_all(d.root, keep_highlight=False)
T = [BM.ptext(p) for p in d.paragraphs()]
pickle.dump(T, open(H / "base_T.pkl", "wb"))
keys = sys.argv[1:]
for k, t in enumerate(T):
    if any(x in t for x in keys):
        print(f"==== base para [{k}]"); print(t)
        for l, f, n in BM.EDITS:
            if f and f in t: print("   EDIT:", l, "| find:", f[:120], "| new:", (n or "")[:300])
        for l, a, n in BM.REPLACE_PARAGRAPHS:
            if t.startswith(a): print("   REPLACE_PARAGRAPH:", l)
