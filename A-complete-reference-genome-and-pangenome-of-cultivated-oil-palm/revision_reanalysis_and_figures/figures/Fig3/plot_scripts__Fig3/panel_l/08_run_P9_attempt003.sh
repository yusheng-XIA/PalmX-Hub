#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/run_config.sh"
assert_run_root
[[ -s "${RUN_CHECKPOINTS}/P9_TABLES_AND_FIGURE_ATTEMPT002.PASS" ]] || { echo "[ERROR] missing accepted P9 attempt002" >&2; exit 73; }
[[ ! -e "${RUN_RESULTS}/P9_final_attempt003_favorable_rescue" ]] || { echo "[ERROR] attempt003 exists" >&2; exit 74; }
touch "${RUN_LOGS}/P9_render_attempt003.out" "${RUN_LOGS}/P9_render_attempt003.err"
"${PYTHON}" -B "${RUN_SCRIPTS}/07_render_favorable_rescue.py" \
  --base-targets "${RUN_RESULTS}/P9_final_attempt002/actionable_targets_with_plot_codes.tsv" \
  --candidates "${RUN_RESULTS}/comprehensive_attempt002/all_trait_relevant_candidates.tsv" \
  --ase "${ASE}" --bins "${BINS}" \
  --output-dir "${RUN_RESULTS}/P9_final_attempt003_favorable_rescue" \
  > "${RUN_LOGS}/P9_render_attempt003.out" 2> "${RUN_LOGS}/P9_render_attempt003.err"
touch "${RUN_PROVENANCE}/P9_scripts.attempt003.sha256" "${RUN_PROVENANCE}/accepted_P9_outputs.attempt003.sha256" "${RUN_CHECKPOINTS}/P9_FAVORABLE_RESCUE_ATTEMPT003.PASS"
sha256sum "${RUN_SCRIPTS}/07_render_favorable_rescue.py" "${RUN_SCRIPTS}/08_run_P9_attempt003.sh" > "${RUN_PROVENANCE}/P9_scripts.attempt003.sha256"
find "${RUN_RESULTS}/P9_final_attempt003_favorable_rescue" -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum > "${RUN_PROVENANCE}/accepted_P9_outputs.attempt003.sha256"
date '+%F %T %z' > "${RUN_CHECKPOINTS}/P9_FAVORABLE_RESCUE_ATTEMPT003.PASS"
echo "P9_ATTEMPT003_PASS"
