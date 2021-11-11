#########
From PyPi
#########

.. _installation:

By far the easiest way to get access to `Dune` is by installing into a
virtual environment using `pip`. In some folder e.g. `Dune` setup a virtual
environment and activate it

.. code-block:: bash

  python3 -m venv dune-env  
  source dune-env/bin/activate

Then download and build `Dune` (note that this takes some time to have
coffee/tea at hand):

.. code-block:: bash

  pip install dune-fem

To test that everything works you can download all the scripts described
in this tutorial and try them out

.. code-block:: bash

  python -m dune.fem
  cd fem_tutorial
  python concepts.py

All examples are available as Python scripts or IPython notebooks - for
this `jupyter` is needed

.. code-block:: bash

  pip install jupyter
  jupyter notebook

This has been tested with different Linux distributions, MacOS, and using
the Windows Linux Subsystem. 

.. note::
   the current installation requires MPI to be available.

###########
From Source
###########

.. note::
   We strongly encourage the use of a python virtual environment and the
   following instructions are written assuming that a virtual environment is
   activated.

************
Requirements
************

The following dependencies are needed for Dune-Fem python binding:

* At least C++17 compatible C++ compiler (e.g. g++ 7 or later)
* python (3.7 or later)

* Required Dune modules (release 2.9 or later)

  * dune-common (https://gitlab.dune-project.org/core/dune-common.git)
  * dune-geometry (https://gitlab.dune-project.org/core/dune-geometry.git)
  * dune-grid (https://gitlab.dune-project.org/core/dune-grid.git)
  * dune-fem (https://gitlab.dune-project.org/dune-fem/dune-fem.git)

* Recommended Dune modules (releases 2.9 or later)

  * dune-istl (https://gitlab.dune-project.org/core/dune-istl.git)
  * dune-localfunctions (https://gitlab.dune-project.org/core/dune-localfunctions.git)
  * dune-alugrid  (https://gitlab.dune-project.org/extensions/dune-alugrid.git)

**********************************
Building the Required Dune Modules
**********************************

Read the instructions on how to `build Dune with Python support`_ which also
links to general instructions on how to `build Dune modules`_.

.. _build Dune modules: https://dune-project.org/doc/installation
.. _build Dune with Python support: https://dune-project.org/doc/pythonbindings

Test your completed installation by opening a Python terminal and running

.. code:: python

   import math
   from dune.grid import structuredGrid
   from dune.fem.function import gridFunction
   grid = structuredGrid([0,0],[1,1],[10,10])
   @gridFunction(grid,name="test",order=2)
   def f(x):
      return math.sin(x.two_norm*2*math.pi)
   f.plot()

If you have everything set up correctly (and have `matplotlib`) you should
get a colored figure and are hopefully ready to go...

.. note::
   The first time you construct an object of a specific realization of one
   of the Dune interfaces (e.g. here a structured grid),
   the just in time compiler needs to be invoked. This can take quite some
   time - especially for grid realizations. This needs to be done only once
   so rerunning the above code a second time (even using other parameters
   in the `structuredGrid` function) should execute almost instantaniously.

***************
Troubleshooting
***************

* The compiler version needs to be 7 or later. This can be checked in terminal with ::

  $ g++ --version

  If your version is out of date, you will need to upgrade your system to use Dune

* It is possible that the python version may be an issue. The scripts
  require `python3` including the development package being installed.
  If during the Dune installation you get the error

  .. code-block:: none

    fatal error: pyconfig.h: No such file or directory

  This can probably be fixed by installing additional python3.5 libraries with e.g. ::

  $ sudo apt-get install libpython3-dev

* One other problem is that a default version of Open MPI may already be installed. This will lead to errors where Dune appears to be looking in the wrong directory for Open MPI (e.g. usr/lib/openmpi instead of the home directory where the script installs it). This can be solved by running ::

  $ make uninstall

  in the original MPI install directory, followed by removing the folder. It will then be necessary to reinstall Open MPI and Dune. It may also be necessary to direct mpi4py to the new MPI installation. It is possible to check whether this is a problem by running python and trying out 

  .. code-block:: python

    from mpi4py import MPI

  If it comes up with an error, this can be fixed by installing mpi4py manually using the following commands ::

  $ git clone https://bitbucket.org/mpi4py/mpi4py.git
  $ cd mpi4py
  $ python setup.py build --mpicc=/path/to/openmpi/bin/mpicc
  $ python setup.py install --user
