#!/usr/bin/env python3
"""Source Data (Fig. 5h,i) and Supplementary Data 24 for the W = 4, P = 15 donor path.

Outputs (this directory):
  sd/Fig5h_chromosome_summary.tsv, sd/Fig5h_window_path.tsv, sd/Fig5h_favourable_loci.tsv,
  sd/Fig5i_selected_chr01B.tsv, sd/Fig5i_load_matrix.tsv      full-sheet replacements for fix_sd_add.py
  (sd/SF10a_constrained_designs.tsv, sd/SF10b_weight_sweep.tsv are written by SF10/make_sf10.py)
  Source_Data_Fig5_W4_candidate.xlsx                            review copy of the split Fig. 5 workbook
  table_cell_changes_5hi_w4.tsv                                 SD24 cell changes (apply_cell_changes.py format)
  fix_sd_add_W4.sh                                              the fix_sd_add.py call (not run here)
Reads the master workbooks read-only; every current W = 0 value is checked against the W = 0 tables first.
"""
import csv, re, shutil
from pathlib import Path
import openpyxl
import pandas as pd

H = Path(__file__).resolve().parent; FIX = H.parent; S = FIX.parent
SD = S / "deliver/Source_Data/Source_Data.xlsx"
SPLIT = S / "deliver/Source_Data_split/Source_Data_Fig5.xlsx"
ST = S / "deliver/Supplementary_Tables/Supplementary_Tables.xlsx"
OUT = H / "sd"; OUT.mkdir(exist_ok=True)


def disp(d):
    m = re.fullmatch(r"(dura|pisifera|nrly|bk)_hap([12])", d)
    return f"{dict(dura='TK', pisifera='NS', nrly='Nigerian', bk='TN')[m.group(1)]}-Hap{m.group(2)}" if m else d


pw = pd.read_csv(H / "path_W4.tsv", sep="\t")
bc = pd.read_csv(H / "by_chrom_W4.tsv", sep="\t")
cp = pd.read_csv(H / "capture_W4_path.tsv", sep="\t")
pc = pd.read_csv(H / "per_chrom_capture_W4.tsv", sep="\t").set_index("Chrom")
w0 = pd.read_csv(FIX / "enh_C/data/ideal_loadonly_path.tsv", sep="\t")
pw["Donor"] = pw.Donor_ID.map(disp)

wb = openpyxl.load_workbook(SD, read_only=True)
def sheet(n): return [list(r) for r in wb[n].iter_rows(values_only=True)]

# ---- sanity: the master holds the W = 0 path
cur = sheet("Fig5h_window_path")
assert len(cur) == 3487 and all(r[4] == disp(d) for r, d in zip(cur[1:], w0.Donor_ID)), "master is not the W = 0 path"

def write(name, header, rows):
    with open(OUT / f"{name}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header)
        for r in rows: w.writerow(["" if v is None or (isinstance(v, float) and pd.isna(v)) else v for v in r])

# 1. chromosome summary (same columns as the current sheet + GWAS weight)
hdr = ["Chrom", "Window_Count", "Residual_DSV", "Residual_DSNP", "Residual_Total_Load", "Breakpoint_Count", "Segment_Count",
       "Donor_Count", "Path_Objective", "Breakpoint_Penalty_P", "GWAS_Weight_W", "Fav_Captured", "Fav_Total",
       "Fav_Capture_Rate_pct", "Fav_No_Donor_Call"]
rows = []
for r in bc.itertuples():
    p = pc.loc[r.Chrom]
    rate = None if p.Fav_Total == 0 else round(100 * p.Captured_W4 / p.Fav_Total, 2)
    rows.append([r.Chrom, r.Window_Count, r.Residual_DSV, r.Residual_DSNP, r.Residual_Total_Load, r.Breakpoint_Count,
                 r.Segment_Count, r.Donor_Count, r.Objective, 15, 4, int(p.Captured_W4), int(p.Fav_Total), rate,
                 int(p.Unknown_W4_no_donor_call)])
write("Fig5h_chromosome_summary", hdr, rows)

# 2. window path
hdr = ["Chrom", "Window_Index", "Window_Start_0based", "Window_End_0based", "Donor_ID", "DSV_Count", "DSNP_Count",
       "Total_Load", "Fav_Loci_Carried", "Donor_Switch_From_Previous", "Cumulative_Chrom_Objective",
       "Cumulative_Genome_Objective", "Breakpoint_Penalty_P", "GWAS_Weight_W"]
write("Fig5h_window_path", hdr, pw[["Chrom", "Window_Index", "Window_Start_0based", "Window_End_0based", "Donor", "DSV_Count",
      "DSNP_Count", "Total_Load", "Fav_Loci_Carried", "Donor_Switch_From_Previous", "Cumulative_Chrom_Objective",
      "Cumulative_Genome_Objective", "Breakpoint_Penalty_P", "GWAS_Weight"]].values.tolist())

# 3. favourable loci (row order as the current sheet)
cur = sheet("Fig5h_favourable_loci"); cpi = cp.set_index("SV")
assert cur[0] == ["SV", "Chrom", "Pos", "Target", "Selected_Donor", "Captured", "Donor_Genotype_Call"]
rows = []
for r in cur[1:]:
    c = cpi.loc[r[0]]; assert int(c.Pos) == r[2] and c.Target == r[3]
    rows.append([r[0], r[1], r[2], r[3], disp(c.W4_Selected_Donor), int(c.Captured_W4),
                 "No_call" if c.Capture_Status_W4.startswith("Unknown") else "Called"])
write("Fig5h_favourable_loci", cur[0], rows)
fav_rows = rows

# 4. chr01B selected windows
cur = sheet("Fig5i_selected_chr01B"); z = pw[pw.Chrom == "chr01B"]
write("Fig5i_selected_chr01B", cur[0], z[["Chrom", "Window_Index", "Window_Start_0based", "Window_End_0based", "Donor",
                                           "DSV_Count", "DSNP_Count", "Total_Load"]].values.tolist())

# 5. 5i load matrix: only Selected_on_path changes
cur = sheet("Fig5i_load_matrix"); assert cur[0][-1] == "Selected_on_path"
sel = {(int(r.Window_Start_0based), r.Donor) for r in z.itertuples()}
rows = [r[:-1] + [int((r[3], r[1]) in sel)] for r in cur[1:]]
assert sum(r[-1] for r in rows) == 354
write("Fig5i_load_matrix", cur[0], rows)
n_sel_changed = sum(1 for a, b in zip(cur[1:], rows) if a[-1] != b[-1])

# ---- Contents descriptions and the fix_sd_add call
DESC = {
    "Fig5h_chromosome_summary": "Per-chromosome summary of the donor path (W = 4, P = 15), with favourable-locus capture",
    "Fig5h_window_path": "Window-level donor path (W = 4, P = 15; 3,486 windows); objective = load - 4 x favourable loci carried + 15 x switches",
    "Fig5h_favourable_loci": "284 SV-GWAS favourable loci and their capture by the W = 4, P = 15 path",
    "Fig5i_selected_chr01B": "chr01B windows selected by the W = 4, P = 15 path",
    "Fig5i_load_matrix": "Load of 35 donor haplotypes in 354 chr01B windows; colour = ln(1 + dSV + dSNP); Selected_on_path for the W = 4, P = 15 path",
}
with open(H / "fix_sd_add_W4.sh", "w") as fh:
    fh.write("#!/bin/bash\n# Replace the Fig. 5h,i sheets with the W = 4 path and add the Supplementary Fig. 10 sheets.\n"
             "# Run from the scratchpad root (same convention as final_rebuild.sh). Not run by the candidate build.\nset -euo pipefail\n"
             "D=fix/fig5hi_w4/sd\npython3 deliver_build/fix_sd_add.py deliver/Source_Data/Source_Data.xlsx \\\n")
    for n, d in DESC.items():
        fh.write(f'  "{n}=$D/{n}.tsv|replace={n}|desc={d}" \\\n')
    fh.write('  "SF10a_constrained_designs=$D/SF10a_constrained_designs.tsv|after=SF9e_QQ|fig=Supplementary Fig. 10|panel=a|'
             'desc=Feasibility-constrained minimum-burden donor designs (African35; at most k breakpoints per chromosome, or at most m donors genome-wide with at most 2 breakpoints per chromosome) and the Fig. 5h path: residual load and reduction relative to the best single donor" \\\n')
    fh.write('  "SF10b_weight_sweep=$D/SF10b_weight_sweep.tsv|after=SF10a_constrained_designs|fig=Supplementary Fig. 10|panel=b|'
             'desc=GWAS-weight sweep (W = 0, 1, 2, 4, 8; P = 15): residual load, breakpoints and favourable loci captured"\n')
os_ = (H / "fix_sd_add_W4.sh"); os_.chmod(0o755)

# ---- review copy of the split workbook
cand = H / "Source_Data_Fig5_W4_candidate.xlsx"; shutil.copy(SPLIT, cand)
wbc = openpyxl.load_workbook(cand)
from openpyxl.styles import Font
for n in DESC:
    i = wbc.sheetnames.index(n); del wbc[n]; ws = wbc.create_sheet(n, i)
    for k, row in enumerate(csv.reader(open(OUT / f"{n}.tsv"), delimiter="\t")):
        ws.append([None if v == "" else (int(v) if re.fullmatch(r"-?\d+", v) else (float(v) if re.fullmatch(r"-?\d+\.\d+", v) else v)) for v in row])
        if k == 0:
            for c in ws[1]: c.font = Font(bold=True)
    cont = wbc["Contents"]
    for r in range(1, cont.max_row + 1):
        if cont.cell(r, 1).value == n: cont.cell(r, 4).value = DESC[n]
wbc.save(cand)

# ---- Supplementary Data 24 cell changes (apply_cell_changes.py format)
st = openpyxl.load_workbook(ST, read_only=True)["Supplementary Table 24"]
T = [list(r) for r in st.iter_rows(values_only=True)]
def cell(ref):
    m = re.fullmatch(r"([A-Z]+)(\d+)", ref); col = openpyxl.utils.column_index_from_string(m.group(1))
    r = T[int(m.group(2)) - 1]; return r[col - 1] if col - 1 < len(r) else None
changes = []
def add(ref, new, need="required", note=""):
    old = cell(ref)
    try: same = old is not None and new is not None and abs(float(old) - float(new)) < 1e-9
    except (TypeError, ValueError): same = False
    if not same and str(old if old is not None else "") != str(new if new is not None else ""):
        changes.append(["Supplementary_Data_24.xlsx", "Supplementary Data 24", ref, "" if old is None else old,
                        "" if new is None else new, need, note])
tot = bc[bc.Chrom == "TOTAL"].iloc[0]; P = pc.loc["TOTAL"]
load, bp, seg, obj, cap = int(tot.Residual_Total_Load), int(tot.Breakpoint_Count), int(tot.Segment_Count), int(tot.Objective), int(P.Captured_W4)
red = 58857 - load; pct = round(100 * red / 58857, 2); rate = round(100 * cap / 284, 2); unk = int(P.Unknown_W4_no_donor_call)
fixed = [
    ("A1", "Supplementary Data 24. Donor path (W = 4, P = 15), candidate deleterious load and favourable-locus capture"),
    ("A2", "Donor-path result (W = 4, P = 15)"),
    ("A4", "African35 donor path (W = 4, P = 15)"),
    ("F4", 4),
    ("A7", "Inputs and exact results (W = 4, P = 15)"),
    ("B17", 4), ("D17", "Reward of 4 per favourable locus carried by the selected donor in the window (genotype-called donors); objective = load - 4 x favourable loci carried + 15 x breakpoints"),
    ("B20", int(tot.Residual_DSNP)), ("B21", load), ("B22", bp), ("B23", seg), ("B24", int(tot.Donor_Count)),
    ("B27", red), ("D27", f"58,857 - {load:,}"), ("B28", pct), ("D28", f"{red:,} / 58,857 × 100"),
    ("B29", cap), ("B31", rate), ("D31", f"{cap} / 284 × 100"),
    # row 32 is a merged note row (A32:N32): one text cell only
    ("A32", f"Note: {unk} favourable loci whose selected donor haplotype has no genotype call are counted as not captured."),
    ("A34", "Chromosome-level summary (W = 4, P = 15)"),
    ("A54", "Window-level donor path, W = 4, P = 15 (3,486 rows)"),
    ("A3544", "Locus-level favourable-allele capture on the W = 4, P = 15 path (284 rows)"),
]
for ref, new in fixed: add(ref, new)
# rows 36-51 per chromosome, 52 total (columns B..M as the sheet: Window_Count .. Fav_Capture_Rate_pct)
for i, r in enumerate(bc.itertuples()):
    rr = 36 + i; assert cell(f"A{rr}") == r.Chrom
    p = pc.loc[r.Chrom]
    vals = dict(C=r.Residual_DSV, D=r.Residual_DSNP, E=r.Residual_Total_Load, F=r.Breakpoint_Count, G=r.Segment_Count,
                H=r.Donor_Count, I=r.Objective, K=int(p.Captured_W4), L=int(p.Fav_Total),
                M=(None if p.Fav_Total == 0 else round(100 * p.Captured_W4 / p.Fav_Total, 2)))
    if r.Chrom == "TOTAL": vals["M"] = rate
    for c, v in vals.items(): add(f"{c}{rr}", None if v is None or pd.isna(v) else (round(float(v), 2) if c == "M" else int(v)))
# window rows 56..3541: Donor_ID (E), DSV (F), DSNP (G), Total (H), switch (I), cum chrom (J), cum genome (K)
assert T[54][:12] == ["Chrom", "Window_Index", "Window_Start_0based", "Window_End_0based", "Donor_ID", "DSV_Count", "DSNP_Count",
                      "Total_Load", "Donor_Switch_From_Previous", "Cumulative_Chrom_Objective", "Cumulative_Genome_Objective", "Breakpoint_Penalty_P"]
for i, r in enumerate(pw.itertuples()):
    rr = 56 + i; assert cell(f"A{rr}") == r.Chrom and cell(f"B{rr}") == r.Window_Index
    for c, v in dict(E=r.Donor, F=r.DSV_Count, G=r.DSNP_Count, H=r.Total_Load, I=r.Donor_Switch_From_Previous,
                     J=r.Cumulative_Chrom_Objective, K=r.Cumulative_Genome_Objective).items():
        add(f"{c}{rr}", v)
# locus rows 3546..3829
for i, r in enumerate(fav_rows):
    rr = 3546 + i; assert cell(f"A{rr}") == r[0]
    add(f"E{rr}", r[4]); add(f"F{rr}", r[5]); add(f"G{rr}", r[6])
with open(H / "table_cell_changes_5hi_w4.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["file", "sheet", "cell", "old_value", "new_value", "necessity", "note"]); w.writerows(changes)
print("sd sheets written; 5i Selected_on_path changed cells:", n_sel_changed)
print("SD24 cell changes:", len(changes))
print(f"TOTAL load {load} bp {bp} seg {seg} obj {obj} capture {cap}/284 = {rate}% reduction {pct}% unknown {unk}")
