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

import numpy
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
#from dune.grid import structuredGrid as leafGridView
from dune.alugrid import aluCubeGrid as leafGridView
#from dune.alugrid import aluSimplexGrid as leafGridView
from dune.common import FieldVector
from dune.grid import reader
from dune.fem import parameter
from dune.fem.space import dgonb as dgSpace
from dune.fem.space import finiteVolume as fvSpace
from dune.fem.operator import molGalerkin as molGalerkin
#from dune.fem.scheme import galerkin as scheme
from dune.fem.scheme import molGalerkin as scheme
from dune.fem.function import uflFunction, integrate
from dune.ufl import Constant
from ufl import TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal
from ufl import dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional, sin, cos
from ufl import as_vector, as_tensor, avg, jump, dS, CellVolume, FacetArea, atan
from dune.femdg.rk import ssp3

Stepper = ssp3(4)

class MolOp:
    def __init__(self,form,space,cfl=0.45):
        self._op = molGalerkin( form )
        self.cfl = cfl
        self.space = space
        self.localTimeStepEstimate = [dt.value/self.cfl]
    def __call__(self,u,v):
        self._op(u,v)
        self.localTimeStepEstimate = [dt.value/self.cfl]
        return
    def stepTime(self,t0, dt0):
        # set time to model time if available
        # since time is not in the form we don't need this here
        pass
        #if hasattr(self._op.model,"time"):
        #    print("Model has time")
        #elif hasattr(self._op.model,"t"):
        #    print("Model has t")

    def applyLimiter(self, u):
        pass

parameter.append({"fem.verboserank": 0})

domain = (reader.dgf, "periodic2.dgf")
#gridView = leafGridView([0, 0], [1, 1], [20, 20])
gridView = leafGridView(domain, dimgrid=2 )
space    = dgSpace(gridView, order=2, storage="fem")
#space  = fvSpace(gridView)

u    = TrialFunction(space)
v    = TestFunction(space)
n    = FacetNormal(space)
he   = avg( CellVolume(space) ) / FacetArea(space)
hbnd = CellVolume(space) / FacetArea(space)
x    = SpatialCoordinate(space)
time = Constant(0., "time")
dt   = Constant(0.001,"dt")


center  = as_vector([ 0.5,0.5 ])
x0 = x[0] - center[0]
x1 = x[1] - center[1]

ux = -4.0*x1
uy =  4.0*x0

# diffusion factor
epsilon = 0.001

# transport direction and upwind flux
b =  as_vector([ux,uy])

def u0(x,t):
    sig2 = 0.004
    sig2PlusDt4 = sig2+(4.0*eps*t)
    xq = ( x0*cos(4.0*t) + x1*sin(4.0*t)) + 0.25
    yq = (-x0*sin(4.0*t) + x1*cos(4.0*t))
    return (sig2/ (sig2PlusDt4) ) * exp (-( xq*xq + yq*yq ) / sig2PlusDt4 )

# upwind (same as LLF in this case)
hatb = (dot(b, n) + abs(dot(b, n)))/2.0

# boundary values (for left/right boundary)
dD   = 1 # conditional((1+x[0])*(1-x[0])<1e-10,1,0)
g    = 0 # conditional(x[0]<0,atan(10*x[1]),0)

# diffusion factor
eps = Constant(epsilon,"eps")
# characteristic function for left/right boundary
dD   = 1
# conditional((1+x[0])*(1-x[0])<1e-10,1,0)
# penalty parameter
beta = Constant( 10*space.order**2 if space.order > 0 else 1,"beta")
exact = uflFunction(gridView, name="exact", order=3, ufl=u0(x,time))

aInternal     = dot(eps*grad(u) - b*u, grad(v)) * dx
diffSkeleton  = eps*dot(avg(grad(v)),jump(u,n))*dS \
               +eps*dot(avg(grad(u)),jump(v,n))*dS \
               -eps*beta/he*jump(u)*jump(v)*dS
diffSkeleton -= eps*beta/hbnd*(u-exact)*v*dD*ds -\
                eps*dot(grad(exact),n)*v*dD*ds
advSkeleton   = jump(hatb*u)*jump(v)*dS
advSkeleton  += hatb*u*v*ds

# minus since we are solving d_t u = L[u]
# which leads in the simplest form to
# unew = uold + dt * L[uold]
# However, in dune-fem we assume that L is on the left hand side
# therefore we implement -L[u]
form = -(aInternal + advSkeleton)
if abs(epsilon) > 1e-12:
    form += diffSkeleton

op = MolOp(form,space,cfl=0.25)
stepper = Stepper(op, cfl=op.cfl)

uh = space.interpolate( exact, name="solution")
uh_n = uh.copy()

# time derivative
explicit = True
endTime = 0.75
if not explicit:
    form -= inner(u - uh_n, v)/dt * dx
    scheme = scheme(form==0, solver="gmres",
                    parameters={"newton.linear.preconditioning.method":"jacobi",
                                "newton.verbose":False,
                                "newton.linear.verbose":False,
                                "newton.linear.tolerance": 1e-12,
                                "newton.tolerance":1e-8
                               })

# <markdowncell>
# <codecell>

error0 = 0.
eoc = 0
for i in range(2):
    time.value = 0.
    uh.interpolate( exact )
    uh.plot()
    gridView.writeVTK("adv",number=0,celldata=[uh],pointdata=[uh,exact])
    while time.value < endTime:
        if explicit:
            stepper( uh, dt.value )
        else:
            uh_n.assign( uh )
            scheme.solve(target=uh)
        time.value += dt.value

    gridView.writeVTK("adv",number=1,celldata=[uh],pointdata=[uh,exact])
    error1= numpy.sqrt( integrate(gridView,dot(uh-exact,uh-exact),order=6))
    uh.plot()
    exact.plot()

    if i > 0:
        eoc = [ numpy.log(error1/error0) / numpy.log(0.5) ]
    print(i,error1, eoc)
    error0 = error1
    gridView.hierarchicalGrid.globalRefine(1)
    dt.value *= 0.5
