#!/bin/bash
# finalize.sh N TYPE BEFORE_PDF : docx 2400-px PNG, before/after comparison, deterministic checks + critic
S=${WORK_DIR}
R=$S/fix/redraw2/ED; N=$1; F=Extended_Data_Fig_0$N
python3 -c "import sys; sys.path.insert(0,'$R/common'); import figtools; figtools.docx_png('$R/out/$F.png','$R/out/docx_2400/$F.png')"
python3 $R/common/figtools.py compare "$3" $R/out/$F.pdf ED$N
python3 $R/common/figtools.py check $R/out/$F.pdf $R/out/$F.png ED$N $2
