import sys, csv, openpyxl
from pathlib import Path
S = Path(__file__).resolve().parents[2]
wb = openpyxl.load_workbook(S / "deliver/Source_Data/Source_Data.xlsx", read_only=True)
rows = list(wb["SF14c_per_accession"].iter_rows(values_only=True)); h = rows[0]
new = list(csv.DictReader(open(S / "fix/sf_plan/sf14/source_data/SF14c_per_accession.tsv"), delimiter="\t"))
inst = {r[0]: r for r in rows[1:]}
assert len(inst) == len(new) == 308
nd = ng = 0
for x in new:
    r = inst[x["Sample_ID"]]
    for i, c in enumerate(h):
        v, w = x[c], r[i]
        if c == "K4_group":
            ng += (str(w) != v)
            continue
        try: nd += abs(float(v) - float(w)) > 1e-9
        except (TypeError, ValueError): nd += (str(v) != str(w))
print(sys.argv[1], "other-column differences:", nd, "| K4_group cells differing from new:", ng)
assert nd == 0
if sys.argv[1] == "post": assert ng == 0
