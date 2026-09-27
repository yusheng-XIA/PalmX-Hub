#!${DATA_DIR}/miniconda3/bin/python
"""Centralize figures, validate expected outputs, and write a concise status report."""
from pathlib import Path
import shutil
import pandas as pd

RUN = Path(__file__).resolve().parents[1]


def main():
    tasks = pd.read_csv(RUN / "manifests/sv_tasks.tsv", sep="\t")
    snp_rows = []; missing = []
    for r in tasks.itertuples(index=False):
        sdir = RUN / "snp/results" / r.category / r.trait; vdir = RUN / "sv/results" / r.category / r.trait
        sp = sdir / "summary_row.tsv"
        if sp.is_file(): snp_rows.append(pd.read_csv(sp, sep="\t"))
        else: missing.append(str(sp))
        for modality, source in [("SNP", sdir), ("SV", vdir)]:
            target = RUN / "figures" / modality / r.category / r.trait; target.mkdir(parents=True, exist_ok=True)
            for name in ["manhattan.png", "manhattan.pdf", "qq.png", "qq.pdf"]:
                src = source / name
                if src.is_file() and src.stat().st_size: shutil.copy2(src, target / name)
                else: missing.append(str(src))
    if snp_rows: pd.concat(snp_rows, ignore_index=True).to_csv(RUN / "snp/summary.tsv", sep="\t", index=False)
    sv = pd.read_csv(RUN / "sv/summary.tsv", sep="\t") if (RUN / "sv/summary.tsv").is_file() else pd.DataFrame()
    audit = pd.read_csv(RUN / "manifests/phenotype_filter_audit.tsv", sep="\t")
    known = pd.read_csv(RUN / "tables/literature_supported_known_gene_overlaps.tsv", sep="\t")
    with (RUN / "reports/FINAL_STATUS.txt").open("w") as fh:
        fh.write(f"Eligible traits after zero and 2.5% two-tail exclusion: {(audit.status == 'eligible').sum()}\n")
        fh.write(f"Skipped traits: {(audit.status != 'eligible').sum()}\n")
        fh.write(f"SNP summaries: {len(snp_rows)}/{len(tasks)}\nSV summaries: {len(sv)}/{len(tasks)}\n")
        fh.write(f"Literature-supported candidate-gene rows: {len(known)}\nMissing required files: {len(missing)}\n")
        if missing: fh.write("\n".join(missing) + "\n")
    if missing: raise RuntimeError(f"final validation failed: {len(missing)} files missing")
    print(f"complete traits={len(tasks)} snp={len(snp_rows)} sv={len(sv)} known_gene_rows={len(known)}")


if __name__ == "__main__": main()
