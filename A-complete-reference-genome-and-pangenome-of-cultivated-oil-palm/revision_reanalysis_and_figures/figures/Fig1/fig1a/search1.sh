B=${ANALYSIS_DIR}/22_answer_reviews
cd $B
echo "== name search"
find . -maxdepth 6 \( -iname '*fao*' -o -iname '*descals*' -o -iname '*oil_palm*cover*' -o -iname '*global*distrib*' -o -iname '*fig*1a*' -o -iname '*figure1*' -o -iname '*.tif' -o -iname '*harvest*' -o -iname '*vegoil*' -o -iname '*contribution*' \) 2>/dev/null | grep -v '/\.' | head -80
echo "== content search"
grep -rIl --include='*.R' --include='*.py' --include='*.md' --include='*.txt' --include='*.tsv' --include='*.csv' --include='*.sh' -e 'FAOSTAT' -e '39\.5%' -e 'Descals' -e 'oil-palm cover' -e 'Oil-palm cover' . 2>/dev/null | head -50
