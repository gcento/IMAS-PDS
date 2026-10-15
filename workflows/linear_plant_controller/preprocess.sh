#!/bin/bash
# The initial state is the same as plasmaless_controller's: NICE inverse equilibrium and
# pf_active of a completed metis_nice_inverse_from_dina run for this shot, in the window
# source.t_min/t_max of this case. Run that workflow's preprocess, which reads this
# case's own workflow_settings.ymmsl and scenario_settings.ymmsl (see its header for the
# PDS_REPO/SCENARIOS_REPO/SHOT/CASE_DIR contract).
exec "$PDS_REPO/workflows/plasmaless_controller/preprocess.sh" "$@"
