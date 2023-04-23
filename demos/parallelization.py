# %% [markdown]
# # Parallelization
# .. index:: Parallelization
#
# Parallelization is available using either distributed memory based on **MPI**
# or multithreading using **OpenMP**.
#
# ## OpenMP
# .. index:: Parallelization; OpenMP
#
# It is straightforward to enable some multithreading support. Note that
# this will speedup assembly and evaluation of the spatial operators but
# not in general the linear algebra backend so that the overall speedup of
# the code might not be as high as expected. Since we rely mostly on
# external packages for the linear algebra, speedup in this step will
# depend on the multithreading support available in the chosen linear
# algebra backend - see the discussion on how to [switch between linear solver backends](solvers_nb.ipynb)
# to, for example, use the thread parallel solvers from scipy.
#
# By default only a single thread is used. To enable multithreading simply add

# %%
from dune.fem import threading
print("Using",threading.use,"threads")
threading.use = 4 # use 4 threads
print("Using",threading.use,"threads")

# %% [markdown]
# At startup a maximum number of threads is selected based on the hardware
# concurrency. This number can be changed by setting the environment
# variable `DUNE_NUM_THREADS` which sets both the maximum and set the number of
# threads to use. To get this number or set the number of threads to use to
# the maximum use

# %%
print("Maximum number of threads available:",threading.max)
threading.useMax()
print("Using",threading.use,"threads")

# %% [markdown]

# ## MPI
# .. index:: Parallelization; MPI
#
# .. todo:: add some remarks on `dune.fem.loadBalance` and possibly some information on loadbalancing strategies available in `dune-alugrid`?
#
#
# It is straightforward to use **MPI** for parallelization. It requires a
# parallel grid, in fact most of the DUNE grids work in parallel except `albertaGrid` or `polyGrid`.
# Most iterative solvers in DUNE work for parallel runs. Some of the
# preconditioning methods also work in parallel, a complete list is found at the [bottom of the solver discussion](solvers_nb.ipynb).
#
# Running a parallel job can be
# done by using `mpirun`
# ```
# mpirun -np 4 python script.py
# ```
# in order to use `4` MPI processes. Example scripts that run in parallel are,
# for example, the [Re-entrant Corner Problem](laplace-adaptive_nb.ipynb).
#
# .. note:: On a cluster where multiple parallel jobs are run simultaneously, it's
# advisable to use one separate cache per job. This can be easily
# done by copying an existing cache with pre-compiled modules. See this [SLURM batch script](slurmbatchscript.rst)
# for an example on how to do this.
#
# ### Load balancing
#
# When running distributed memory jobs load balancing is an issue. There are two
# ways to ensure a balanced work load.

# %%

from dune.alugrid import aluCubeGrid as leafGridView
from dune.grid import cartesianDomain
domain = cartesianDomain([0,0], [1,1], [40,40], overlap=1)
gridView = leafGridView( domain )
# print possible load balancing methods
help( leafGridView )
