Tutorial for the Python Bindings of DUNE-FEM Module
===================================================

*Authors: Andreas Dedner, Martin Nolte, and Robert Kl&ouml;fkorn*

This repository for the **dune-fem** modules provides the Python scripts that form the
[tutorial][tutlink] for the [dune-fem][femlink] discretization
module. Here is an example:

```Python
# from laplace-readme.py

import numpy as np
from ufl import *
import dune.ufl
from dune.ufl import DirichletBC
from dune.grid import cartesianDomain
from dune.alugrid import aluConformGrid
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import uflFunction

domain = cartesianDomain([0, 0], [1, 1], [16, 16])
gridView = aluConformGrid( domain, dimgrid=2)

# set up a space and a conforming finite element scheme and solve the PDE
space  = lagrange(gridView, dimRange=1, order=1)
uh = space.function(name="u_h")

# set up a diffusion reaction model using UFL
u = TrialFunction(space)
v = TestFunction(space)
x = SpatialCoordinate(space)

# provide an exact solution that will be used to add suitable forcing and Dirichlet b.c.
exact = as_vector( [cos(2.*pi*x[0])*cos(2.*pi*x[1])] )
# make 'exact' into a grid function for output and uh into an UFL coefficient for error computation
exact_gf = uflFunction( gridView, "exact", order=5, ufl=exact)

# bilinear form
a = (inner(grad(u), grad(v)) + inner(u,v)) * dx
# forcing term (in this case very simple)
b = (inner(grad(exact_gf), grad(v)) + inner(exact_gf,v)) * dx

# create a scheme (operator and solver)
dbc = DirichletBC(space, exact)
scheme = galerkin([a==b, dbc])

# compute solution
info = scheme.solve( target=uh )

# now define a grid function representing the pointwise error
l2error_gf = uflFunction( gridView, "error", 5, as_vector([dot((exact-uh),(exact-uh))]) )

error = np.sqrt( l2error_gf.integrate() )
print("size:",gridView.size(0),"L2-error:",error)
gridView.writeVTK("laplace", pointdata=[ uh, l2error_gf, exact_gf ])
```

See the file COPYING for full copying permissions.

If you find this [tutorial][tutlink] helpful for getting your own projects up and
running please cite

Title: Python Bindings for the DUNE-FEM module
*Authors: Andreas Dedner, Martin Nolte, and Robert Klöfkorn*
Publisher: Zenodoo, 2020
DOI 10.5281/zenodo.3706994

[![DOI](https://zenodo.org/badge/246695081.svg)](https://zenodo.org/badge/latestdoi/246695081)

Dependencies
------------

dune-fempy depends on the following DUNE modules:
- dune-common 2.9+
- dune-geometry 2.9+
- dune-grid 2.9+
- dune-python 2.9+
- dune-fem 2.9+

In addition to the dependencies of these DUNE modules, the following software
packages are required:

- a C++17 compatible C++ compiler
- Python 3.6+

We strongly recommend installing the following Python packages to make full
use of dune-fem's Python bindings:

- numpy
- mpi4py
- ufl (2022.2.0+)

[femlink]: https://gitlab.dune-project.org/dune-fem/dune-fem
[tutlink]: https://dune-project.org/sphinx/content/sphinx/dune-fem/
