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
