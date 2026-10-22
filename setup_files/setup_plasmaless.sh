# this file expects to be run from the 'local_installs' folder
set -euo pipefail # stop if anything doesn't work

# Linear plant models, MATLAB (actor muscle3/muscle_linear_model_actor.m, plant code
# plant/, models linear_models/, machine description data/md_dd4/). Nothing to build: the
# install is the checkout. MATLAB (with the Control System Toolbox), PCS (for its muscle3
# venv) and IMAS-MATLAB are loaded at run time by the `linear_model` program in
# workflows/lib/local_programs.ymmsl, not here.
#
# BRANCH_PLASMALESS defaults to the IMAS_Muscle3 branch of the fork, which holds the
# MUSCLE3 actor with the linear models (main does not). pds_setup.sh passes it too.
# A tag or commit may be passed as 2nd argument instead to pin a version.
PLASMALESS_URL=${1:-"https://github.com/gcento/plasma-tokamak-circuits.git"}
BRANCH_PLASMALESS=${2:-"IMAS_Muscle3"}

if [[ ! -d "plasmaless-tokamak-circuits/.git" ]]; then
  git clone $PLASMALESS_URL plasmaless-tokamak-circuits
fi
cd plasmaless-tokamak-circuits
git fetch --quiet origin
git checkout $BRANCH_PLASMALESS
echo "  plasmaless commit: $(git rev-parse --short HEAD)"
cd ..
