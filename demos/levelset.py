# %% [markdown]
# .. index:: Equations; Levelset Advection
#
# # Levelset Advection
# %%
try:
    import dune.femdg
except ImportError:
    print("This example needs 'dune.femdg'. Install dune-fem-dg with pip or by source!")
    import sys
    sys.exit(0)

import numpy
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from ufl import *
from dune.grid import cartesianDomain
from dune.alugrid import aluSimplexGrid as simplexGrid
from dune.alugrid import aluCubeGrid
from dune.fem.space import dglagrangelobatto, finiteVolume
from dune.fem.function import gridFunction
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
from dune.fem.utility import lineSample
from dune.ufl import cell, Constant

dune.fem.threading.use = 4

# dimension
dim = 2

# interface thinckness
eps = 0.01

# %% [markdown]
# .. index:: Methods; Discontinuous Galerkin - Spectral Element Method (DG-SEM)
#
# Using a DG-Lagrange space with Gauss-Lobatto points as interpolation nodes we obntain
# a discontinuous Galerkin - Spectral Element method (DG-SEM).

def getSpace( gv, order, dgsem=True ):
    return dglagrangelobatto(gridView, dimRange=1, order=order ) if dgsem else \
           dgonb( gridView, dimRange=1, order=order )

# %% [markdown]
# Basic model for advection
# %%
class Model:
    dimRange = 1

    # min and max values of U
    lowerBound = [0.0]*dimRange
    upperBound = [1.0]*dimRange

    endTime = 8

    _gamma = None # gamma # diffusion coefficient
    _eps = eps # interface thickness

    def setVelocity( v ):
        try:
            x = SpatialCoordinate(cell(dim))
            v(0, x)
            Model._v = v
        except (ValueError, TypeError):
            Model._v = lambda x,t: v

    def F_c(t,x,U):
        v = Model.velocity(t,x,U)
        return as_matrix( [[ *(v*U[0]) ]])

    # interface re-initialization anti diffusion
    if _gamma is not None:
        def coeff(U, DU=None ):
            norm_DU = 1.0 if DU is None else sqrt(inner(DU,DU) + 1e-8)
            return Model._gamma*(Model._eps) #  - U[0]*(1-U[0])/norm_DU)
        def F_v(t,x,U, DU):
            print("Diffusion enabled!")
            return Model.coeff( U, DU ) * DU
        def maxDiffusion(t,x,U):
            return Model.coeff( U )

    # max wave speed of advection
    def maxWaveSpeed(t,x,U,n):
        v = Model.velocity(t,x,U)
        return abs(dot(v,n))

    # velocity method
    def velocity(t,x,U):
        if Model._v is not None:
            return as_vector(Model._v(t,x))
        else:
            T = Model.endTime
            u = -2*sin(pi*x[0])*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[1])*cos(pi*t/T)
            v =  2*cos(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(pi*t/T)
            return as_vector([u,v])

    # variable should be between lowerBound and upperBound
    def physical(t,x,U):
        return ( conditional( U[0] < Model.upperBound[0] + 1e-6, 1.0, 0.0 ) *
                 conditional( U[0] > Model.lowerBound[0] - 0.001, 1.0, 0.0 ) )

    # jump of the variable
    def jump(t,x,U,V):
        return U-V

    # indicator for adaptation
    class Indicator:
        def eta(t,x,U):
            return U[0]
        def F(t,x,U,DU):
            q = Model.F_c(t,x,U)[0]
            return q
        def S(t,x,U,DU):
            return 0

# %% [markdown]
# Method to evolve the solution in time
# %%
def evolve(space, u_h, limiter=None, T=8):
    lu_h = u_h.localFunction()
    #@gridFunction(space.gridView,name="phi",order=space.order)
    operator = femDGOperator(Model, space, limiter=limiter)
    stepper  = femdgStepper(order=space.order+1 if space.order>0 else 2, operator=operator)
    # does nothing here
    operator.applyLimiter(u_h)

    t = 0
    saveStep = 0.1
    #fig = pyplot.figure(figsize=(30,10))
    u_h.plot() # gridLines="white", figure=(fig, 131), colorbar=False, clim=[0,1])
    c = 1
    while t <= T:
        operator.setTime(t)
        t += stepper(u_h)
        if t > saveStep:
            print(t)
            saveStep += 0.1
            #u_h.plot() # gridLines="white", figure=(fig, 131+c), colorbar=False, clim=[0,1])
            c += 1

    u_h.plot() # gridLines="white", figure=(fig, 131+c), colorbar=False, clim=[0,1])
    res = numpy.zeros((2,space.gridView.size(0)))
    for i,e in enumerate(space.gridView.elements):
        x = e.geometry.center
        res[0][i] = x.two_norm
        res[1][i] = lu_u(e,e.geometry.toLocal(x))
    return res

# %% [markdown]
# We start with a triangular grid
# %%

domain = cartesianDomain([0]*dim,[1]*dim, [20]*dim, periodic=[False]*dim, overlap=1)
gridView = aluCubeGrid( domain, dimgrid=2 )
#gridView.plot()

# fix the order to use and storage structure for results
res = {}
order = 2
space = getSpace( gridView, order )

# %% [markdown]
# A time reversed rotating vortex
# %%
x = SpatialCoordinate(cell(dim))

T = 8

# velocity for the time reversed vortex
def v(t,x):
    u = -2*sin(pi*x[0])*sin(pi*x[0])*cos(pi*x[1])*sin(pi*x[1])*cos(pi*t/T)
    v =  2*cos(pi*x[0])*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[1])*cos(pi*t/T)
    return [u,v]

Model.setVelocity( v )

# center of circle
center = [0.5, 0.75]
radius = 0.15

Model.epsilon = Constant(0.01, "eps")

# circle as initial condition
diff = sqrt(((x[0] - center[0])**2 + (x[1] - center[1])**2))

phi = radius - diff - 2*Model.epsilon

psi = 0.5*tanh( phi/(2*Model.epsilon) + 1 ) + 0.5
initial = as_vector([psi])


# %% [markdown]
# Initial conditions for a radial Riemann problem
#
# First solved on a triangle grid and minmod limiter

# %%
u_h   = space.interpolate( initial, name="solution")
res["simplex (minmod)"] = evolve(space,u_h)
