# %% [markdown]
# .. index:: Boundary; Dirichlet Conditions
#
# # More General Boundary Conditions
# So far we only used natural boundary conditions. Here we discuss how to
# set Dirichlet boundary conditions and use different conditions for
# different components of the solution.
#
# To fix Dirichlet boundary conditions $u=g$ on part of the boundary
# $\Gamma\subset\partial\Omega$ the central class is
# `dune.ufl.DirichletBC` which takes three arguments:
# the discrete function space for $u$, the function $g$ given by a UFL
# expression, and a description of $\Gamma$. There are different ways to
# do this. If it is omitted or `None` the condition is applied to the whole
# domain, a integer $s>0$ can be provided which can be set to describe a
# part of the boundary during grid construction as described at the end of
# this section. Finally a UFL condition can be used, i.e., `x[0]<0`.
#
# For vector valued functions $u$ the value function $g$ can be a UFL
# vector or a list. In the later case a component of `None` can be used to
# describe components which are not to be constrained by the boundary
# condition.

# %%
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import ticker
from dune.fem.plotting import plotComponents

import numpy, io
from dune.grid import structuredGrid as leafGridView
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.function import gridFunction, cppFunction, integrate, uflFunction
from dune.ufl import DirichletBC, Constant
from ufl import TestFunction, TrialFunction, SpatialCoordinate,\
                dx, ds, grad, inner, sin

gridView = leafGridView([0, 0], [1, 1], [4, 4])
vecSpace = solutionSpace(gridView, dimRange=2, order=2)
x = SpatialCoordinate(vecSpace)
vec = vecSpace.interpolate([0,0], name='u_h')
uVec,vVec = TrialFunction(vecSpace), TestFunction(vecSpace)
a  = ( inner(grad(uVec), grad(vVec)) + inner(uVec,vVec) ) * dx
f  = ( uVec[0]*(1-uVec[1])*vVec[0] + uVec[1]*(1-uVec[0])*vVec[1] ) * dx
# Add nonlinear Neuman boundary conditions for all non Dirichlet boundaries
# for second quantity:
f  = f + uVec[0]*uVec[0] * vVec[1] * ds
# Define dirichlet Boundary conditions for first component on all boundaries
bc = DirichletBC(vecSpace,[sin(4*(x[0]+x[1])),None])
vecScheme = solutionScheme( [a == f, bc],
        parameters={"newton.linear.tolerance": 1e-9} )

vecScheme.solve(target=vec)

plotComponents(vec, gridLines=None, level=2,
               colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})

# %% [markdown]
# To prescribe $u_2=0$ at the bottom boundary is also straightforward - the
# Neuman boundary flux added to the right hand side `f` previously will
# still be used on the top and vertical boundaries:

# %%
bcBottom = DirichletBC(vecSpace,[sin(4*(x[0]+x[1])),1],x[1]<1e-10)
vecScheme = solutionScheme( [a == f, bc, bcBottom],
        parameters={"newton.linear.tolerance": 1e-9} )
vecScheme.solve(target=vec)
plotComponents(vec, gridLines=None, level=2,
               colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})

# %% [markdown]
# We can also use general grid functions (including discrete functions)
# for the boundary conditions in the same way we could use grid functions
# anywhere within the UFL forms. So the following code leads to the same
# results as above:

# %%
test = vec.copy()
test.clear()
bc = DirichletBC(vecSpace,[vec[0],None])
@gridFunction(gridView ,name="bnd",order=2)
def bnd(x):
    return [numpy.sin(4*(x[0]+x[1])),1]
bcBottom = DirichletBC(vecSpace,bnd,x[1]<1e-10)
vecScheme = solutionScheme( [a == f, bc, bcBottom],
        parameters={"newton.linear.tolerance": 1e-9} )
vecScheme.solve(target=test)
plotComponents(test, gridLines=None, level=2,
               colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})
assert sum([ abs(td-vd) for td,vd in zip(test.dofVector, vec.dofVector)] ) /\
       len(vec.dofVector) < 1e-7

test.clear()
code="""
#include <cmath>
template <class GridView>
auto bnd() {
  return [](const auto& en,const auto& xLocal) -> auto {
    auto x = en.geometry().global(xLocal);
    return sin(4.*(x[0]+x[1]));
  };
}
"""
bndCpp = cppFunction(gridView, name="bndCpp", order=2,
                     fctName="bnd",includes=io.StringIO(code))
bc = DirichletBC(vecSpace,[bndCpp,None])
vecScheme = solutionScheme( [a == f, bc, bcBottom],
                            parameters={"newton.linear.tolerance": 1e-9} )
vecScheme.solve(target=test)
plotComponents(test, gridLines=None, level=2,
               colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})
assert sum([ abs(td-vd) for td,vd in zip(test.dofVector, vec.dofVector)] ) /\
       len(vec.dofVector) < 1e-7

# %% [markdown]
# We provide a few method, to apply or otherwise access information about
# the dirichlet constraints from an existing  scheme.

# The following shows how to set all degrees of freedom on the Dirichlet
# boundary to zero.

# %%
dc = vecScheme.dirichletBlocks
for i,block in enumerate(dc):
    for j,b in enumerate(block):
        if b > 0:
            vec.dofVector[i*len(block)+j] = 0
plotComponents(vec, gridLines=None, level=2,
               colorbar={"orientation":"horizontal", "ticks":ticker.MaxNLocator(nbins=4)})

# %% [markdown]
# Up until now we have used ufl boolean expressions to
# define the subset of the boundary where a Dirichlet condition should be
# applied. In many cases different parts of the boundary are already marked
# with ids (integers) during grid construction. These ids can also be used
# to define where boundary conditions are to be applied.
# At the moment to use boundary ids the grid has to be constructed using
# the Dune Grid Format (dgf) reader
#
# todo: add tutorial link or describe here and link to dgf documentation).
# todo: describe what the boundary ids for a Cartesian domain are.
#
# These ids can be used as third argument in the `DirichletBC` constructor
# and also in the boundary measure `ufl.ds` as shown in the following
# examples solving a Laplace problem with a combination of Neumann, Robin, and
# Dirichlet boundary conditions.

# %% [markdown]
# In this example the boundary conditions are
#
# - left and right boundary (id=1 and id=2):
#   Robin boundary conditions $\nabla u\cdot n + u = -10$ for $y\in[0.4,0.6]$ and
#   natural boundary conditions ($\nabla u\cdot n=0$) otherwise.
# - bottom boundary (id=3):
#   Neumann conditions $\nabla u\cdot n=2$
# - top boundary (id=4):
#   Dirichlet conditions $u=1$ for $x\in [0.3,0.6]$ and
#   natural boundary conditions ($\nabla u\cdot n=0$) otherwise.

# %%
from dune.ufl import BoundaryId
from ufl import ds, grad, dot, inner, conditional, eq, And, FacetNormal, Not

gridView.hierarchicalGrid.globalRefine(3)
space = solutionSpace(gridView, order=2)
x,n = SpatialCoordinate(space), FacetNormal(space)
u,v = TrialFunction(space), TestFunction(space)
a  = inner(grad(u), grad(v)) * dx

fbnd  = (u-10 * v *                    # middle of left and right boundary
         conditional(abs(x[1]-0.5)<0.1,1,0) ) * ds( (1,2) )
fbnd += 2 * v * ds(3)                  # bottom boundary
dbc   = DirichletBC(space, 1,          # middle of top boundary
           And( eq(BoundaryId(space),4), abs(x[0]-0.5)<0.2 ))

scheme   = solutionScheme( [a == fbnd, dbc] )
solution = space.interpolate(0, name='u_h')
scheme.solve(target=solution)
solution.plot()

# %% [markdown]
# The final example has a slightly more complex setup at the boundary:
#
# - left and right (id=1 and id=2):
#   for $y\in[0.4,0.6]$ we set Neumann conditions $\nabla u\cdot n=-10$
# - rest of boundary: 1

# %%
fbnd  = (-10 * v *                     # middle of left and right boundary
         conditional(abs(x[1]-0.5)<0.1,1,0) ) * ds( (1,2) )
dbc   = DirichletBC(space, 1,        # not middle of left and right boundary
           Not( And(BoundaryId(space)<=2,abs(x[1]-0.5)<0.1) ) )

scheme   = solutionScheme( [a == fbnd, dbc] )
solution = space.interpolate(0, name='u_h')
scheme.solve(target=solution)
solution.plot()
