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
cd doc
# - make
# make -k ${parallel_opts} concepts_nb.ipynb
python concepts.py
cd ..
${DUNECONTROL} --current make doc
