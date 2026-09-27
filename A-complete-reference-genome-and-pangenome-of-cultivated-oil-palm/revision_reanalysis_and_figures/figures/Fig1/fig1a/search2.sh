B=${ANALYSIS_DIR}/22_answer_reviews/00_ms
cd $B
which pdftotext
for f in 03_V3/01_figure1/Fig1-7.28_1.pdf "03_V3/01_figure1/Fig1-7.28(1).pdf" 05_MS/new_revision/Final_figures/Figure1-A4_1.5fold.pdf 05_MS/new_revision/Final_figures/figure1-5.pdf 05_MS/0918_revision/Figure1a_oil_palm_cultivation_points.pdf 05_MS/0918_revision/投稿材料整合_待核对_20260923/01_正文图/Figure1.pdf; do
 echo "=== $f"; pdftotext -l 1 "$f" - 2>/dev/null | grep -i -n "%\|FAO\|harvest\|oil output\|cover\|Soy\|Rape\|Sunflower\|Contribution\|major" | head -40
done
echo "=== ref mapping"
grep -i "FAOSTAT\|8\.5%\|8\.6%\|39\.5\|39%" 04_最终正文/01_revision_ms/参考文献逐条核验与主张映射_20260805.tsv | head -20
