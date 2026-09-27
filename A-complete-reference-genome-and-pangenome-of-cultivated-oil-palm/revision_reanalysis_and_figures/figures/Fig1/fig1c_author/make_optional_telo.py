"""Optional supplement (not a replacement): Fig1c_telomeres with per-end telomere call under the final rule.
Data rows = our current Fig1c_telomeres (identical to author sheet 7); call from fix/telomere/telo_ends.tsv."""
import csv, openpyxl
from pathlib import Path
S=Path(__file__).resolve().parents[2]
calls={}; win={}
for r in csv.DictReader(open(S/'fix/telomere/telo_ends.tsv'),delimiter='\t'):
    if r['label'] in ('EG11','FL-Hap2') and r['seq'].startswith('chr') and not r['seq'].endswith('B'):
        calls[(r['label'],r['seq'],r['end'])]=int(r['sd2_call_ge3']); win[(r['label'],r['seq'],r['end'])]=r['max_win']
ws=openpyxl.load_workbook(S/'deliver/Source_Data_split/Source_Data_Fig1.xlsx',read_only=True)['Fig1c_telomeres']
rows=list(ws.iter_rows(values_only=True))
out=[[rows[0][0]],list(rows[1])+['Telomere-positive (final rule)','Max repeats in a terminal 10-kb window']]
n={}
for r in rows[2:]:
    if r[0] in ('EG11','FL-Hap2'):
        k=(r[0],r[1],'L' if r[4]=='left' else 'R'); c=calls[k]; n.setdefault(r[0],0); n[r[0]]+=c
        out.append(list(r)+['Yes' if c else 'No',int(win[k])])
out.append([])
out.append([rows[-1][0]+' Telomere-positive: at least one 10-kb window within the terminal 50 kb contains three or more TTTAGGG/CCCTAAA copies (Methods; Supplementary Data 2). Positive ends: EG11 %d/32, FL-Hap2 %d/32.'%(n['EG11'],n['FL-Hap2'])])
with open(Path(__file__).with_name('Fig1c_telomeres_OPTIONAL_with_calls.tsv'),'w',newline='') as f:
    csv.writer(f,delimiter='\t',lineterminator='\n').writerows(out)
print(n)
