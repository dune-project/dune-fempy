# <markdowncell>
# # Advection-Diffusion: Discontinuous Galerkin Method with Upwinding
# So far we have been using Lagrange spaces of different order to solve our
# PDE. In the following we show how to use Discontinuous Galerkin method to
# solve an advection dominated advection-diffusion probllem:
# \begin{align*}
# -\varepsilon\triangle u + b\cdot\nabla u &= f
# \end{align*}
# with Dirichlet boundary conditions. Here $\varepsilon$ is a small
# constant and $b$ a given vector.
# <codecell>

import numpy, math
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.grid import structuredGrid as leafGridView
from dune.common import FieldVector
from dune.fem import parameter
from dune.fem.space import dglegendre as dgSpace
from dune.fem.space import finiteVolume as fvSpace
from dune.fem.scheme import molGalerkin as molScheme
from dune.fem.scheme import galerkin as galerkinScheme
from dune.fem.function import uflFunction, integrate
from dune.ufl import Constant
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal
from ufl import dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional, sin, cos
from ufl import as_vector, avg, jump, dS, CellVolume, FacetArea, atan

parameter.append({"fem.verboserank": 0})

gridView = leafGridView([-1, -1], [1, 1], [20, 20])
space    = dgSpace(gridView, order=0, storage="istl")
#space  = fvSpace(gridView)

u    = TrialFunction(space)
v    = TestFunction(space)
n    = FacetNormal(space)
he   = avg( CellVolume(space) ) / FacetArea(space)
hbnd = CellVolume(space) / FacetArea(space)
x    = SpatialCoordinate(space)
time = Constant(0., "time")
dt = Constant(0.001,"dt")

# diffusion factor
eps = Constant(0.0,"eps")
# transport direction and upwind flux
b    = as_vector([1,1])

#hatb = avg(u) + 1/2 * jump(u)
hatb = (dot(b, n) + abs(dot(b, n)))/2.0
flux = hatb*u
#hatb = conditional( dot(b, n) < 0, dot(b, n), 0 )
# boundary values (for left/right boundary)
dD   = 1 # conditional((1+x[0])*(1-x[0])<1e-10,1,0)
g    = 0 # conditional(x[0]<0,atan(10*x[1]),0)
eps = None
center  = as_vector([ 0.0,0.0 ])

aInternal     = dot(- b*u, grad(v)) * dx
advSkeleton   = jump(flux)*jump(v)*dS
advSkeleton  += flux*v*ds

xt = x - time*b
exact = conditional( dot(xt,xt) < 0.05, 1, 0 )
exactFct =  uflFunction(gridView, name="exact", order=1, ufl=exact)

uh = space.interpolate( exactFct, name="solution")
uh_n = uh.copy()
uh_old = uh.copy()

# time derivative
explicit = True
endTime = 0.5
form  = dt * (aInternal + advSkeleton)
if explicit:
    form -= inner(u, v ) * dx
else:
    form += inner(u - uh_n, v) * dx

scheme = molScheme(form==0, solver=("suitesparse","umfpack"),
                   parameters={"newton.linear.preconditioning.method":"jacobi",
                               "newton.verbose":True,
                               "newton.linear.verbose":True,
                               "newton.tolerance":1e-4
                              })

# <markdowncell>
# <codecell>

error0 = 0.
eoc = 0
for i in range(2):
    time.value = 0.
    uh.interpolate( exact )
    uh.plot()
    while time.value < endTime:
        uh_n.assign( uh )
        if explicit:
            scheme(uh_n, uh )
        else:
            scheme.solve(target=uh)
        time.value += dt.value

    uh.plot()
    error1= math.sqrt( integrate(gridView,dot(uh-exactFct,uh-exactFct),order=2))
    if i > 0:
        eoc = [ math.log(error1/error0) / math.log(0.5) ]
    print(i,error1, eoc)
    error0 = error1
    gridView.hierarchicalGrid.globalRefine(1)
    dt.value *= 0.5
