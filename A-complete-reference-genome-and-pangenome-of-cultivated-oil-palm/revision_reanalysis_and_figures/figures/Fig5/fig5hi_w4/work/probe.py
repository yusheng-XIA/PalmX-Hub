import re, sys
from pathlib import Path
S = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(S / "deliver_build"))
import build_main_text as B
from docx_accept import accept_all
base = S / "v1639/HB-01_MS_FinalFigures_HighRes_Legends_Tracked.docx"
d = B.Doc(base); accept_all(d.root, keep_highlight=False)
texts = [B.ptext(p) for p in d.paragraphs()]
keys = ["W = 4", "346", "362 donor", "60.04", "23,510", "F_{GWAS}", "favourable-locus score", "reward term", "Haplotype-informed donor prioritization", "exact dynamic programming", "Several limitations define"]
for i, t in enumerate(texts):
    if any(k in t for k in keys):
        print(f"\n=== BASE [{i}] {t}")
        for lab, f, n in B.EDITS:
            if f and f in t: print(f"   EDIT {lab!r}: find@{t.find(f)}-{t.find(f)+len(f)} :: {f[:120]!r}")
        for lab, anc, n in B.REPLACE_PARAGRAPHS:
            if t.startswith(anc): print("   REPLACE_PARAGRAPH", lab)
