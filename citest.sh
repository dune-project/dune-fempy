#! /bin/bash
set -e
set -u
set -- --opts="/duneci/dune.opts" "${@}"
DUNECONTROL=dunecontrol

pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple -r reqCI.txt
cd /duneci/modules
./dune-common/bin/setup-dunepy.py
cd $CI_PROJECT_DIR

parallel_opts=
if [[ -v DUNECI_PARALLEL ]]; then
  echo "Parallel run with ${DUNECI_PARALLEL} processes"
  parallel_opts="-j${DUNECI_PARALLEL}"
fi
export OMPI_MCA_rmaps_base_oversubscribe=1
export OMPI_MCA_mpi_yield_when_idle=1
export OMPI_MCA_btl_base_warn_component_unused=0
set -x

duneci-standard-test
cd doc
# make -k ${parallel_opts} concepts_nb.ipynb
cd ..
${DUNECONTROL} --current make doc
