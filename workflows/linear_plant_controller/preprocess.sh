#!/bin/bash
# The initial state is built as plasmaless_controller's: NICE inverse equilibrium and
# pf_active of preprocess.source_uri (this workflow's default: the inverse_convergence run
# of the shot, settings.ymmsl), in the window source.t_min/t_max of this case. Run that
# workflow's preprocess, which reads both from this case's own workflow_settings.ymmsl and
# scenario_settings.ymmsl (see its header for the PDS_REPO/SCENARIOS_REPO/SHOT/CASE_DIR
# contract).
exec "$PDS_REPO/workflows/plasmaless_controller/preprocess.sh" "$@"
