#!/bin/bash
#
#SBATCH -A <PROJECT-ID>
#SBATCH -t 00:10:00
#SBATCH -J <JOB-NAME>
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err
#
# Number of nodes the job should use
#SBATCH -N 1
# Number of cores per node
#SBATCH --tasks-per-node=4

###
# Warning! This template script assumes:
# * your Dune installation is in `$DUNE_HOME`
# * you installed Dune using the `dune-fem-dg/scripts/build-dune-fem-dg.sh` script
# * you want to create the `$RUNDIR` directory to store job data
# * This should also work fine with other installations, i.e. pip install but
# * might require some adjustments.
###
DUNE_HOME="$HOME/dune"
RUNDIR="`pwd`/jobs/$SLURM_JOB_ID"

: '
Optional comment describing the purpose of the run.
'

# To create a separate venv may not be necessary but it should not hurt either
python3 -m venv dune-env
# activate virtual environment
source dune-env/bin/activate

# Use local compilation cache.
# To use the global cache - comment the following lines.

# copy global cache to reuse compiled modules
cp -r "$DUNE_HOME/cache" "$RUNDIR"
# set environment variable to new cache directory
export DUNE_PY_DIR="$RUNDIR/cache"
# force re-configuration of cache to adjust paths
rm "$DUNE_PY_DIR/dune-py/.noconfigure"

echo "JOB STARTED"

# Example, replace this part with the test case you intended to run.
mpirun \
python $DUNE_HOME/testcase.py

echo "JOB FINISHED"
