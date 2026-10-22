#!/bin/bash
# Builds the initial state (the F_INIT equilibrium + pf_active of magnetic_controller and
# of the plant) of linear_plant_controller once, at case creation, instead of inside every
# run. It replaces the former in-workflow chain source (nice_out) -> waveform_editor,
# making the same calls (preprocess_initial_state.py) with this directory's
# waveforms.yaml. Run once by bin/pds-create-case, output frozen into $CASE_DIR/preprocess/
# (not rebuilt on every bin/pds-run-case.sbatch). See that script's header for the
# PDS_REPO/SCENARIOS_REPO/SHOT/CASE_DIR contract.
#
# source.source_uri in settings.ymmsl -> ${CASE_DIR}/preprocess/initial_state
#
# Inputs, both read from this case's own workflow_settings.ymmsl + scenario_settings.ymmsl
# (already written into $CASE_DIR by bin/pds-create-case when this runs; the override wins),
# so each is defined in one place only:
#  - preprocess.source_uri: the NICE inverse output the plasma and coil data come from;
#    default in settings.ymmsl: the inverse_convergence run (out_nice) of the shot; a case
#    override can point at another run;
#  - source.t_min/t_max: the window of that output.
set -euo pipefail

WF_DIR="$PDS_REPO/workflows/linear_plant_controller"
OUT="$CASE_DIR/preprocess"
mkdir -p "$OUT"

# Frozen copy of the waveform file, its ${SCENARIOS_REPO}/${SHOT} placeholders resolved
# (machine description from pds-scenarios ${SHOT}/data/in_md), kept for provenance.
envsubst '${PDS_REPO} ${SCENARIOS_REPO} ${SHOT} ${CASE_DIR}' \
  < "$WF_DIR/waveforms.yaml" > "$OUT/waveforms.yaml"

# Same module as the waveform_editor program in workflows/lib/easybuild_programs.ymmsl
# (base_env: clean there, hence the purge). If PDS got loaded as an actual module
# (bin/pds-create-case's own auto-load, or the caller's shell), Lmod adopted PDS_REPO via
# its setenv and purge unwinds that -- reassert it after.
_PDS_REPO="$PDS_REPO"
module purge
module load Waveform-Editor/0.3.2.dev154-intel-2025b-pds
export PDS_REPO="$_PDS_REPO"

# preprocess.source_uri from the stacked case settings (the later file wins), read as
# source_settings() in preprocess_initial_state.py reads the window
SOURCE_URI="$(python - "$CASE_DIR/workflow_settings.ymmsl" "$CASE_DIR/scenario_settings.ymmsl" <<'EOF'
import sys
from pathlib import Path

import yaml

merged = {}
for path in sys.argv[1:]:
    if Path(path).is_file():
        merged.update((yaml.safe_load(Path(path).read_text()) or {}).get("settings") or {})
print(merged.get("preprocess.source_uri", ""))
EOF
)"
if [[ -z "$SOURCE_URI" ]]; then
  echo "preprocess.sh: no preprocess.source_uri in the case settings: set it in the" \
       "workflow's settings.ymmsl or in the case override" >&2
  exit 1
fi
SOURCE_DIR="${SOURCE_URI#"imas:hdf5?path="}"
if [[ "$SOURCE_DIR" != "$SOURCE_URI" && ! -d "$SOURCE_DIR" ]]; then
  echo "preprocess.sh: no NICE inverse output $SOURCE_DIR: run that inverse workflow for" \
       "shot $SHOT first, or point preprocess.source_uri at another run" >&2
  exit 1
fi
echo "preprocess.sh: initial state from $SOURCE_URI"

python -u "$WF_DIR/preprocess_initial_state.py" \
  --source-uri "$SOURCE_URI" \
  --waveforms "$OUT/waveforms.yaml" \
  --out "$OUT/initial_state" \
  "$CASE_DIR/workflow_settings.ymmsl" "$CASE_DIR/scenario_settings.ymmsl"
