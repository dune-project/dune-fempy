.. dune-fem python documentation master file, created by
   Andreas Dedner on Mon Mai 20 2019.

################################
Welcome to the dune-fem tutorial
################################

This module brings python scripting support to `Dune`_.
This version describes the bindings for the development version
(to become 2.10). The bindings serves three purposes:

1. High level program control for solving partial differential equations
   using classes from the `Dune`_ core and from `Dune-Fem`_
   :cite:`Dune-Grid-Paper,Dune-Fem-Paper` with a recent update provided in :cite:`DuneReview`.
   The unified form language `UFL`_ :cite:`UFL-Paper`
   is used to describe the mathematical
   model, all realizations of the `Dune`_ grid interface can be used to
   work with the domain tessellation, and the finite element spaces,
   operator, and solvers provided by `Dune-Fem`_ for the discretizations 
   and solving steps. All of this is available to be used in python
   scripts or through Jupyter notebooks.
2. Rapid prototyping of new methods or other parts of a simulation is easy
   since the interfaces provided are very similar to the `Dune`_ C++
   interface. This makes it easy to transfer a working prototype from
   python (easy to develop) to C++ (high efficiency). Small C++ code
   snippets can be easy called from python using just in time compilation.
3. Rapid prototyping of new implementations of `Dune`_ interfaces. For example
   new implementations of the `Dune`_ grid interface can be easily
   tested. For `Dune-Fem`_ developers, new grid views, discrete function spaces, and
   scheme classes can be added and tested.

.. _Dune: https://www.dune-project.org
.. _Dune-Python: https://www.dune-project.org/modules/dune-python/
.. _Dune-Istl: https://www.dune-project.org/modules/dune-istl/
.. _Dune-Fem: https://www.dune-project.org/modules/dune-fem/
.. _Dune-Vem: https://www.dune-project.org/modules/dune-vem/
.. _Dune-Fem-Dg: https://www.dune-project.org/modules/dune-fem-dg/
.. _UFL: https://bitbucket.org/fenics-project/ufl

.. toctree::
   :maxdepth: 1
   :caption: Overview

   overview

.. toctree::
   :maxdepth: 2
   :caption: Installation

   installation

.. toctree::
   :maxdepth: 3
   :caption: Getting Started

   dune-fempy_nb
   concepts_nb


.. toctree::
   :maxdepth: 2
   :caption: Next Steps

   nextsteps

.. toctree::
   :maxdepth: 2
   :caption: Further Examples

   furtherexamples

.. toctree::
   :maxdepth: 3
   :caption: Further Topics

   furthertopics

.. toctree::
   :maxdepth: 3
   :caption: Extension Modules

   extensionmodules

.. toctree::
   :maxdepth: 1
   :caption: User Projects
   :name: userprojects

   userprojects

.. toctree::
   :maxdepth: 2
   :caption: Information and Resources

   inforesources

.. API reference <api/modules>

.. index::
   see: Mesh construction; Grid construction
   see: Load Balancing; Parallelization

############
Bibliography
############

.. bibliography:: dune-fempy.bib
   :all:
