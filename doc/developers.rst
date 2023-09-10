##############################
Information for C++ Developers
##############################

The `dune` Python bindings rely on just in time (JIT) compilation which
takes place in a special `dune` module `dune-py`. Which each C++ type a
Python module is generated there. One can obtain information about and
manipulate this module using the `dune` main script:

.. code-block::

  python -m dune
  usage: dune [-h] {info,configure,list,remove,fix-dunepy,dunetype,checkbuilddirs}

  remove              Remove generated modules
  info                Print information about dune-py
  list                List all generated modules
  dunetype            Show dune types for given modules

  fix-dunepy          Find inconsistencies in dune-py and try to fix automatically. This will potentially delete all generated modules.
  configure           Tag dune-py to be reconfigured before next use
  checkbuilddirs      (internal use) Check build directories

Probably the most important option here is `remove` which simplifies
removing specific modules from `dune-py`. This is important to debug C++
code because a module can then be recompiled with debugging compiler flags.

.. code-block::

  python -m dune remove
  usage: dune remove [-h] [--beforedate] [modules [modules ...]]

  positional arguments:
    modules       Patterns of modules ("\*.cc" and dune-py path is added to each argument) or "all"

  optional arguments:
    -h, --help    show this help message and exit
    --beforedate  Instead of a pattern provide a date to remove all modules not having been loaded
                  after that date

In addition to providing a full module (without extension and path),
one can provide the module prefix (e.g. `hierarchicalgrid`) and all modules
with this prefix will be removed. If is quite easy to find the module name
from a generated `dune` object in Python by using the `__module__` magic
attribute. The following script outputs

.. code-block::

  from dune.grid import structuredGrid
  view = structuredGrid([0,0],[1,1],[10,10])
  print(view.__module__)

`dune.generated.hierarchicalgrid_966e2a5c8356c5b278ccd3acad180f0a`.
So removing this module can be achieved by calling

.. code-block::

  python -m dune remove hierarchicalgrid_966e2a5c8356c5b278ccd3acad180f0a

There are two options to compile a module that has not been generated yet
(e.g. was removed) with different flags. Either by setting the `CXXFLAGS`
environment variable, e.g.,

.. code-block::

  CXXFLAGS=-g python test.py

Alternatively, we can set the flags from within a Python script

.. code-block::

  import dune.generator as generator
  generator.addToFlags("-DWANT_CACHED_COMM_MANAGER=0",noChecks=True)
  algorithm(...)
  generator.setFlags("-g -Wfatal-errors",noChecks=True)
  algorithm(...)
  generator.reset()

Now one can use `gdb` to debug the C++ code

.. code-block::

  gdb python
  (gdb) r test.py

Other available environment variables are

* `DUNE_LOG_LEVEL`: set to `info` or `debug` to get more verbose output
  during module building.
* `DUNE_SAVE_BUILD`: show full output from compiler during JIT compilation.
  Set to `terminal` to obtain output to console or set to `write` to obtain
  log files in current folder.
* `DUNE_PY_DIR`: location of `dune-py` module
* `DUNE_CMAKE_FLAGS`: cmake flags to use during configuration of a new `dune-py` module
* `DUNEPY_DISABLE_PLOTTING`: disable plotting using the `dune` plotting functions
