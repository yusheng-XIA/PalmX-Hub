#!/bin/bash
set -e
set -u
set -o pipefail

ROOT="${DATA_DIR2}/projects/1-oil_palm/08-Nature_review/7-Fig1/10-syntenyviz-eg11-fl-africa"
PYTHON_BIN="${DATA_DIR2}/anaconda3/bin/python3"
MANIFEST="${ROOT}/config/Input_Manifest.tsv"
PARAMETERS="${ROOT}/config/Parameters.tsv"
SOURCE_PAF="${DATA_DIR2}/projects/1-oil_palm/08-Nature_review/7-Fig1/results/02_synteny/EG11_vs_FL_Africa_hap2/Orientation.paf"
FINAL_STAGE="${ROOT}/results/EG11_FL_synteny_inputs"
ARCHIVE="${ROOT}/deliverables/EG11_FL_synteny_inputs.tar.gz"
ARCHIVE_SHA="${ARCHIVE}.sha256"

printf '[INFO] Read-only preflight started | Host=%s | Time=%s\n' "$(hostname)" "$(date --iso-8601=seconds)"

for executable in "${PYTHON_BIN}" /usr/bin/time tar sha256sum gzip; do
    if ! command -v "${executable}" >/dev/null 2>&1; then
        printf '[ERROR] Missing executable: %s\n' "${executable}" >&2
        exit 1
    fi
done

for required in "${MANIFEST}" "${PARAMETERS}" "${SOURCE_PAF}" "${ROOT}/scripts/build_synteny_inputs.py"; do
    if [[ ! -s "${required}" ]]; then
        printf '[ERROR] Missing or empty input: %s\n' "${required}" >&2
        exit 1
    fi
    stat --printf='[INPUT] %n\t%s bytes\tmtime=%y\n' "${required}"
done

for forbidden in "${FINAL_STAGE}" "${ARCHIVE}" "${ARCHIVE_SHA}"; do
    if [[ -e "${forbidden}" ]]; then
        printf '[ERROR] Refusing overwrite; output already exists: %s\n' "${forbidden}" >&2
        exit 1
    fi
done

"${PYTHON_BIN}" -c 'from pathlib import Path; p=Path("${DATA_DIR2}/projects/1-oil_palm/08-Nature_review/7-Fig1/10-syntenyviz-eg11-fl-africa/scripts/build_synteny_inputs.py"); compile(p.read_text(), str(p), "exec"); print("[PASS] Python syntax")'
"${PYTHON_BIN}" "${ROOT}/scripts/build_synteny_inputs.py" --help >/dev/null

"${PYTHON_BIN}" -c 'import csv,pathlib,sys
p=pathlib.Path(sys.argv[1]); rows=list(csv.DictReader(p.open(newline=""), delimiter="\t"))
assert len(rows)==2 and {r["Sample_ID"] for r in rows}=={"EG11","FL"}
for row in rows:
    sample=row["Sample_ID"]
    for key in ("Genome_FASTA","Annotation_GFF","Chromosome_Map"):
        q=pathlib.Path(row[key]); assert q.is_file() and q.stat().st_size>0, q
        print(f"[INPUT] {sample}:{key}\t{q}\t{q.stat().st_size} bytes")
' "${MANIFEST}"

df -h "${ROOT}"
printf '[PASS] Preflight completed; no source or result file was modified.\n'
