#!/bin/bash
R=${DATA_DIR}/youzong/rawdata
ls $R $R/analysis 2>/dev/null | head -80
echo "## find"
nice -n 19 find $R -maxdepth 6 \( -iname '*.hic' -o -iname '*.cool' -o -iname '*.mcool' -o -iname '*porec*' -o -iname '*pore_c*' -o -iname '*pore-c*' -o -iname '*contact_map*' -o -iname '*readscov*' -o -iname '*gaps_fill*' \) 2>/dev/null | grep -v "/envs/" | head -100
