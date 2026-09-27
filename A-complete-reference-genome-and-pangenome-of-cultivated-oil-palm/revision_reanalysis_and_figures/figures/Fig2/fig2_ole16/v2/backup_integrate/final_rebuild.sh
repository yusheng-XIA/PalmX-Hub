#!/bin/bash
# Final rebuild of every deliverable, with the 2026-09-24/25 restructure (ED/SF renumbering by first citation,
# four figures promoted to ED, ED6f-g deleted, Table 2 -> Supplementary Data 18).  One command, re-runnable.
#
# Contract: all build SOURCES stay in the OLD numbering (legends.py, edits_*.py, si_text.py, Source_Data.xlsx sheet
# names, figure file names).  Renumbering happens only in the outputs:
#   deliver_renum/           renamed figure files, renumbered Source Data (+ split), docx images
#   deliver/01_Main_Text     manuscripts (renumbered at build time)
#   deliver/Supplementary_Information, deliver/Extended_Data_Figures   SI and ED Word/PDF
#   deliver/Supplementary_Data  (Table 2 added to SD18 as a worksheet)
# Any failing step stops the run (set -e); fix/restructure/restructure.py check must pass before assembly.
# (The 2026-09-24 one-off 'fix_sd_add' patch of 4c/ED4b/5c sheets is already in Source_Data.xlsx and is not re-run.)
set -euo pipefail
S=$(cd "$(dirname "$0")/.." && pwd)
cd "$S"
R=$S/fix/restructure; B=$S/deliver_build; D=$S/deliver; N=$S/deliver_renum
BASE=$S/v1639/HB-01_MS_FinalFigures_HighRes_Legends_Tracked.docx
RS="python3 $R/restructure.py"

echo "== 1 plan: trial build in the old numbering, first-citation order -> deliver_build/renumber_map.json"
$RS plan --build "$B" --base "$BASE" --figdir "$D"/Main_Figures_revised | tail -28
echo "== 2 install renumbering hooks (idempotent)"
$RS install --build "$B"
echo "== 3 ED6 (-> ED9) without panels f,g"
$RS ed6 --plot "$S"/fix/redraw2/ED/work/ED6/plot_ED6.py --out "$R"/ed6_no_fg | tail -1
echo "== 4 figure files under the new names -> deliver_renum"
$RS figures --map "$B"/renumber_map.json --deliver "$D" --docx-img "$S"/docx_img --override ED6="$R"/ed6_no_fg --out "$N"
echo "== 5 Supplementary Data (+ Table 2 as a worksheet of Supplementary Data 18)"
python3 "$B"/build_supp_data.py "$D"/Supplementary_Tables/Supplementary_Tables.xlsx "$D"/Supplementary_Data | tail -1
$RS supp-data --map "$B"/renumber_map.json --sd-dir "$D"/Supplementary_Data
echo "== 6 Source Data: rename / drop / renumber, cell-by-cell verification, split"
$RS source-data --map "$B"/renumber_map.json "$D"/Source_Data/Source_Data.xlsx "$N"/Source_Data/Source_Data.xlsx | grep -v " -> "
$RS verify-source-data --map "$B"/renumber_map.json --jobs 6 "$D"/Source_Data/Source_Data.xlsx "$N"/Source_Data/Source_Data.xlsx
rm -f "$N"/Source_Data_split/Source_Data_*.xlsx
python3 "$B"/split_source_data.py "$N"/Source_Data/Source_Data.xlsx "$N"/Source_Data_split | tail -2
echo "== 7 SI / ED documents"
python3 "$B"/build_docx.py "$N"/docx_img "$D"/Supplementary_Information/Supplementary_Information.docx \
  "$D"/Extended_Data_Figures/Extended_Data_Figures_with_legends.docx
soffice --headless --convert-to pdf --outdir "$D"/Supplementary_Information "$D"/Supplementary_Information/Supplementary_Information.docx >/dev/null 2>&1
soffice --headless --convert-to pdf --outdir "$D"/Extended_Data_Figures "$D"/Extended_Data_Figures/Extended_Data_Figures_with_legends.docx >/dev/null 2>&1
echo "== 8 main text (clean, tracked, highlighted, reference)"
python3 "$B"/build_main_text.py "$BASE" "$D"/01_Main_Text "$D"/Main_Figures_revised | tail -3
mkdir -p ms_work/pdf_r4
soffice --headless --convert-to pdf --outdir ms_work/pdf_r4 "$D"/01_Main_Text/Manuscript_clean.docx "$D"/01_Main_Text/Manuscript_tracked_changes.docx >/dev/null 2>&1
echo "== 9 check"
M=$D/01_Main_Text
$RS check --map "$B"/renumber_map.json --build "$B" --manuscript "$M"/Manuscript_clean.docx \
  --variant "$M"/Manuscript_tracked_changes.docx --variant "$M"/Manuscript_highlighted_vs_first_submission.docx \
  --variant "$M"/Reference_all_sentences_changed_vs_first_submission.docx \
  --si "$D"/Supplementary_Information/Supplementary_Information.docx \
  --ed "$D"/Extended_Data_Figures/Extended_Data_Figures_with_legends.docx --split-dir "$N"/Source_Data_split
echo "== 10 records"
P=$(python3 -c "import fitz;print(fitz.open('ms_work/pdf_r4/Manuscript_clean.pdf').page_count, fitz.open('ms_work/pdf_r4/Manuscript_tracked_changes.pdf').page_count)")
python3 "$B"/write_apply_record.py "$M" $P
python3 "$B"/build_changes_doc.py "$D"/正文需修改清单.docx
cp "$B"/decisions.md "$D"/需拍板与确认事项.md; pandoc "$D"/需拍板与确认事项.md -o "$D"/需拍板与确认事项.docx
echo "== 11 assemble"
FIGROOT="$N" SOURCE_DATA_SPLIT="$N"/Source_Data_split bash "$B"/assemble_final.sh | tail -1
echo "== 12 archive"
bash "$B"/build_archive.sh | tail -2
echo DONE
