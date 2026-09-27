#!${DATA_DIR}/miniconda3/bin/python
from pathlib import Path
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
rows = []
for r in tasks.itertuples(index=False):
    p = RUN / "sv/results" / r.category / r.trait / "summary_row.tsv"
    if not p.is_file(): raise FileNotFoundError(p)
    rows.append(pd.read_csv(p, sep="\t"))
out = pd.concat(rows, ignore_index=True)
out.to_csv(RUN / "sv/summary.tsv", sep="\t", index=False)
print(f"sv_summaries={len(out)}")
