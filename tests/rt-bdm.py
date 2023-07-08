# NOTE: there is some issue with failing convergence when using solver=cg -
# it should work...
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
import numpy as np
from dune.grid import cartesianDomain
from dune.alugrid import aluSimplexGrid as aluGridView
from ufl import SpatialCoordinate, CellVolume, TrialFunction, TestFunction,\
                inner, dot, div, grad, dx, as_vector, transpose, Identity, sin, pi
from dune.ufl import Constant, DirichletBC
from dune.fem.plotting import plotPointData as plot
from dune.fem.function import integrate
from dune.fem.scheme import galerkin
from dune.fem.space import raviartThomas  as rtSpace
from dune.fem.space import bdm as bdm
import pprint
from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

order = 1
# Note: structuredGrid fails in precon step!
# grid = structuredGrid([0,0],[3,1],[30,10])
grid = aluGridView( cartesianDomain([0,0],[3,1],[30,10]) )

spcR0 = rtSpace( grid, order=0, storage="numpy" )
spcR1 = rtSpace( grid, order=1, storage="numpy" )
spcB1 = bdm( grid, order=1, storage="numpy" )
spcB2 = bdm( grid, order=2, storage="numpy" )

cell  = spcR0.cell()
x     = SpatialCoordinate(cell)
exact = grad(sin(pi*x[0])*sin(pi*x[1]))
# exact = as_vector([-exact[0],exact[1]])

fig, axs = pyplot.subplots(4, 3, figsize=(40,30))

grid.hierarchicalGrid.globalRefine(1)
# error in u0_y and u1_x
for l in range(2):
    results = []
    # div(u).div(v) = -grad(div(u)).v + div(u) v.n
    # u(x,y) = (pi cos(pi x)sin(pi y) , pi sin(pi x)cos(pi y)
    # div u  = -pi^2 sin(pi x)sin(pi y) - pi^2 sin(pi x)sin(pi y)
    #        = -2pi^2 sin(pi x)sin(pi y)
    for i,spc in enumerate([spcR0,spcR1,spcB1,spcB2]):
        uh = spc.interpolate( exact, name="U" )
        u,v = TrialFunction(spc), TestFunction(spc)
        a  = div(u)*div(v) * dx # + dot( grad(div(exact)),v) * dx
        a += dot(u,v) * dx      # - dot(exact,v) * dx
        dbc = DirichletBC(spc,exact)
        scheme = galerkin([a==0,dbc])
        # scheme.solve(target=uh)
        plot(uh[0],   grid=grid, gridLines=None, figure=(fig, axs[i][0]), colorbar="horizontal")
        plot(uh[1],   grid=grid, gridLines=None, figure=(fig, axs[i][1]), colorbar="horizontal")
        plot(div(uh), grid=grid, gridLines=None, figure=(fig, axs[i][2]), colorbar="horizontal")
        # uh.plot()
        error = uh-exact
        results += [[ np.sqrt(e) for e in integrate(grid,[error**2,div(error)**2], order=5) ]]
    pyplot.show()
    pprint.pprint(results)
    grid.hierarchicalGrid.globalRefine(1)
