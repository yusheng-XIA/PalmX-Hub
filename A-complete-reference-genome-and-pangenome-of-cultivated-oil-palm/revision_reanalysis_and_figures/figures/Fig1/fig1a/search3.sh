for B in ${DATA_DIR}/youzong ${DATA_DIR}; do
echo "##### $B"
timeout 900 find $B -maxdepth 9 \( -iname '*faostat*' -o -iname '*FAOSTAT*' -o -iname 'Production_Crops*' -o -iname '*QCL*' -o -iname '*L2_2019*' -o -iname '*oilpalm_2019*' -o -iname '*oil_palm_2019*' -o -iname '*descals*' -o -iname '*GlobalOilPalm*' -o -iname '*oil*palm*cover*' -o -iname '*op_cover*' -o -iname '*harvested*' -o -iname '*vegetable*oil*' -o -iname '*crop_share*' -o -iname '*fig1a*' -o -iname '*figure1a*' \) -not -path '*/11_map_zone/*' 2>/dev/null | head -100
done
