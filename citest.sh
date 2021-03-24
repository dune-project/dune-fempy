#! /bin/bash
set -e
set -u
set -- --opts="/duneci/dune.opts" "${@}"
DUNECONTROL=dunecontrol

pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple --no-build-isolation wheel
pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple --no-build-isolation --upgrade Pygments
pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple --no-build-isolation -r reqCI.txt
# check including portalocker
# pip install -i https://gitlab.dune-project.org/api/v4/projects/812/packages/pypi/simple --no-build-isolation portalocker
cd /duneci/modules
./dune-common/bin/setup-dunepy.py
cd $CI_PROJECT_DIR

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

ret = 0
duneci-standard-test

cd doc
if [[ $REBUILD == "on" ]]; then
  make clean
  make -ki ${parallel_opts}
  makeResult = `make -n`
  if [ -z "$makeResult" ]; then
    ret = 1
  fi
else
  touch concepts.py
  make concepts_nb.ipynb
  ret = "$?"
fi
cd ..

${DUNECONTROL} --current make doc

if [[ $REBUILD == "on" ]]; then
  tar cvzf doc.tar.gz doc
  curl --header "JOB-TOKEN: $CI_JOB_TOKEN" --upload-file $CI_PROJECT_DIR/doc.tar.gz $REGURL/tutorials/0.0.1/tutorial.tar.gz
fi

exit "$ret"
