#!/bin/bash
# finalize_ed.sh N SRC_PDF SRC_PNG TYPE : copy to ED/, docx png, compare, check
S=${WORK_DIR}
B=$S/fix/beautify; N=$1; F=Extended_Data_Fig_0$N
cp "$2" $B/ED/$F.pdf; cp "$3" $B/ED/$F.png
python3 -c "import sys; sys.path.insert(0,'$B/common'); import figtools; figtools.docx_png('$B/ED/$F.png','$B/ED/docx_2400/$F.png')"
python3 $B/common/figtools.py compare $S/deliver/Extended_Data_Figures/$F.pdf $B/ED/$F.pdf ED$N
python3 $B/common/figtools.py check $B/ED/$F.pdf $B/ED/$F.png ED$N $4
