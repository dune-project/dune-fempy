% dune-fem python documentation master file, created by
% Andreas Dedner on Mon Mai 20 2019.

# Welcome to the dune-fem tutorial (for version 2.10)

This module brings python scripting support to [Dune].
This version describes the bindings for the development version
(to become 2.10). The bindings serves three purposes:

1. High level program control for solving partial differential equations
   using classes from the [Dune] core and from [Dune-Fem]
   {cite}`Dune-Grid-Paper,Dune-Fem-Paper` with a recent update provided in {cite}`DuneReview`.
   The unified form language [UFL] {cite}`UFL-Paper`
   is used to describe the mathematical
   model, all realizations of the [Dune] grid interface can be used to
   work with the domain tessellation, and the finite element spaces,
   operator, and solvers provided by [Dune-Fem] for the discretizations
   and solving steps. All of this is available to be used in python
   scripts or through Jupyter notebooks.
2. Rapid prototyping of new methods or other parts of a simulation is easy
   since the interfaces provided are very similar to the [Dune] C++
   interface. This makes it easy to transfer a working prototype from
   python (easy to develop) to C++ (high efficiency). Small C++ code
   snippets can be easy called from python using just in time compilation.
3. Rapid prototyping of new implementations of [Dune] interfaces. For example
   new implementations of the [Dune] grid interface can be easily
   tested. For [Dune-Fem] developers, new grid views, discrete function spaces, and
   scheme classes can be added and tested.

# Table of Content

```{toctree}
:caption: Overview
:numbered:

overview
```

```{toctree}
:caption: Installation
:maxdepth: 1
:numbered:

installation
```

```{toctree}
:caption: Getting Started
:maxdepth: 3
:numbered:

dune-fempy_nb
concepts_nb
```

```{toctree}
:caption: Next Steps
:maxdepth: 2
:numbered:

nextsteps
```

```{toctree}
:caption: Further Examples
:maxdepth: 2
:numbered:

furtherexamples
```

```{toctree}
:caption: Further Topics
:maxdepth: 2
:numbered:

furthertopics
```

```{toctree}
:caption: Extension Modules
:maxdepth: 1
:numbered:

extensionmodules
```

```{toctree}
:caption: User Projects
:maxdepth: 1
:name: userprojects
:numbered:

userprojects
```

```{toctree}
:caption: Information and Resources
:maxdepth: 2

inforesources
```

```{toctree}
:caption: Additional pages
:hidden: true

gmsh2dgf_nb
slurmbatchscript
scheme_api
todo
```

% API reference <api/modules>

```{index} see: Mesh construction; Grid construction see: Load Balancing; Parallelization see: Models; Equations
```

# Bibliography

```{eval-rst}
.. bibliography:: dune-fempy.bib
   :all:
```

[dune]: https://www.dune-project.org
[dune-fem]: https://www.dune-project.org/modules/dune-fem/
[dune-fem-dg]: https://www.dune-project.org/modules/dune-fem-dg/
[dune-istl]: https://www.dune-project.org/modules/dune-istl/
[dune-python]: https://www.dune-project.org/modules/dune-python/
[dune-vem]: https://www.dune-project.org/modules/dune-vem/
[ufl]: https://bitbucket.org/fenics-project/ufl
