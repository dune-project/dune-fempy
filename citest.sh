#! /bin/bash
set -e
set -u
set -- --opts="/duneci/dune.opts" "${@}"
DUNECONTROL=dunecontrol
parallel_opts=
if [[ -v DUNECI_PARALLEL ]]; then
  echo "Parallel run with ${DUNECI_PARALLEL} processes"
  parallel_opts="-j${DUNECI_PARALLEL}"
fi
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_btl_base_warn_component_unused=0
set -x

cd $CI_PROJECT_DIR
python3 -m venv --system-site-packages fempy
source fempy/bin/activate
pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple -r reqCI.txt
# - /duneci/modules/dune-common/bin/setup-dunepy.py install
cd doc
# - make
# make -k ${parallel_opts} concepts_nb.ipynb
python concepts.py
cd ..
${DUNECONTROL} --current make doc
