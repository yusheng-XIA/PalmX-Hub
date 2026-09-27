#!/bin/bash
D=${ANALYSIS_DIR}/01_seedless_results1/15_hifiasm_wihporec
ls -la $D/cphasing_output/5.plot/ | head -30
ls -la $D | head -60
which python3; python3 -c "import h5py, numpy; print('h5py ok')" 2>&1
for p in ${DATA_DIR}/miniconda3/envs/*/bin/python ${DATA_DIR}/*/envs/*/bin/python; do [ -x $p ] && $p -c "import cooler" 2>/dev/null && echo "COOLER $p"; done 2>/dev/null | head -5
