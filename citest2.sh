#! /bin/bash
# set -e
set -u
set -- --opts="/duneci/dune.opts" "${@}"
DUNECONTROL=dunecontrol

# pip install -i https://gitlab.dune-project.org/api/v4/projects/541/packages/pypi/simple --no-build-isolation wheel
# pip install -i https://gitlab.dune-project.org/api/v4/projects/541/packages/pypi/simple --no-build-isolation --upgrade Pygments
# pip install -i https://gitlab.dune-project.org/api/v4/projects/541/packages/pypi/simple --no-build-isolation --upgrade -r reqCI.txt

# python -m ipykernel install --user --name=dune-env
# jupyter kernelspec list

# cd /duneci/modules
# cd $CI_PROJECT_DIR

echo "*******************************"
pip list
echo "*******************************"

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

export DUNE_LOG_LEVEL=DEBUG
export DUNEPY_DISABLE_PLOTTING=1
set +u
echo "CXXFLAGS=$DUNE_CMAKE_FLAGS PATH=$PWD"
set -u
ret=0
cd demos
python laplace-adaptive.py
if [ $? -ne 0 ] ; then ret=1 ; fi
python solversInternal.py
if [ $? -ne 0 ] ; then ret=1 ; fi

exit "$ret"
