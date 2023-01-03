##############################
Information for C++ Developers
##############################

.. todo:: this section could use more details

Options for recompiling certain modules with different compiler flags

.. code-block::

  import dune.generator as generator
  generator.addToFlags("-DWANT_CACHED_COMM_MANAGER=0",noChecks=True)
  algorithm(...)
  generator.setFlags("-g -Wfatal-errors",noChecks=True)
  algorithm(...)
  generator.reset()

Obtaining information about the JIT `dune` module:

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


Of special interest for developers is the `remove` option:

.. code-block::

  python -m dune remove --help
  usage: dune remove [-h] [--beforedate] [modules [modules ...]]

  positional arguments:
    modules       Patterns of modules (\*.cc and dune-py path is added to each argument) or all
  optional arguments:
    --beforedate  Instead of a pattern provide a date to remove all modules not having been loaded after that date

You can use `gdb` to debug the C++ code

.. code-block::

  gdb python
  (gdb) r script.py

The main issue is to compile the module you would like to debug with `-g`.
You can use 

.. code-block::

  generator.setFlags("-g -Wfatal-errors",noChecks=True)

before the construction of the python object for that module or use
`python -m dune remove modulename` if you know which module is is (or use
`python -m dune remove femspace*` for example to debug a space. Then you
can use

.. code-block::

  CXXFLAGS=-g gdb python
  (gdb) r script.py
