import sys
from pathlib import Path
S = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(S / "deliver_build"))
import build_main_text as B
doc = B.Doc(S / "v1639/HB-01_MS_FinalFigures_HighRes_Legends_Tracked.docx")
B.accept_all(doc.root, keep_highlight=False)
P = [B.ptext(p) for p in doc.paragraphs()]
for idx in map(int, sys.argv[1:]):
    print("=== para", idx)
    for l, f, n in B.EDITS:
        if f in P[idx]:
            s = P[idx].index(f)
            print(f"  [{s}:{s+len(f)}] {l}\n     FIND: {f[:200]}\n     NEW : {str(n)[:300]}")
    for l, a, n in B.REPLACE_PARAGRAPHS:
        if P[idx].startswith(a): print("  REPLACE_PARAGRAPH", l)
