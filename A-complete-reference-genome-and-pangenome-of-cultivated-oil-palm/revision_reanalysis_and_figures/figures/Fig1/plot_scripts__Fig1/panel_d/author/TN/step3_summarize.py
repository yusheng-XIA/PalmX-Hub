# step3_summarize.py
import pandas as pd

def count_sv(bed_file):
    try:
        return len(pd.read_csv(bed_file, sep="\t", header=None))
    except:
        return 0

# 外圈：TE重叠 vs 无TE重叠（颜色：genic/intergenic再细分）
n_te      = count_sv("sv_with_te.bed")
n_no_te   = count_sv("sv_no_te.bed")
n_genic   = count_sv("sv_genic.bed")
n_intergenic = count_sv("sv_intergenic.bed")
total     = count_sv("sv_all.bed")

# 内圈：基因体各区域
n_up2k    = count_sv("sv_upstream2k.bed")
n_exon    = count_sv("sv_exon.bed")
n_intron  = count_sv("sv_intron.bed")
n_down2k  = count_sv("sv_downstream2k.bed")

summary = {
    "total": total,
    "te_overlap": n_te,
    "no_te": n_no_te,
    "genic": n_genic,
    "intergenic": n_intergenic,
    "up2k": n_up2k,
    "exon": n_exon,
    "intron": n_intron,
    "down2k": n_down2k
}

pd.Series(summary).to_csv("sv_summary_for_plot.csv", header=False)
print(pd.Series(summary))
