# %% [markdown]
# # General Concepts (with a time dependent problem)
# In the following we introduce how to construct and use the basic components of
# a finite element method ending the first part with solving a simple Laplace
# problem. The required steps are:
#
# *  constructing a tessellation of the computational domain (**grid**)
# *  working with functions defined over the grid (**grid functions**)
# *  setting up **discrete functions**
# *  defining the mathematical model to solve (**ufl**)
# *  solving the (non linear) system arising from the discretization of the
#    model by the Galerkin method (**scheme**)

# %% [markdown]
# .. index::
#    pair: Grids; Grid Definition
#
# ## Setting up the Grid
#
# .. note:: The goal of a Dune **grid** is to provide an
# interface to describe subdivisions (elements) of the computational domain $\Omega$
# in a generic way. See <cite data-cite="DuneReview"></cite>, Section 3.1 for
# more details.
#
# After some general import statements we start our tutorial by setting up
# a simple Cartesian grid over the domain $\Omega = [0,1]^2$ subdivided into four
# intervals in each dimension. We can then show the grid using the `plot`
# method - note that direct plotting with MatPlotLib is only supported in
# 2D, other options for visualization, e.g., using ParaView will be
# discussed later.
#
# Using the `dune.grid.structuredGrid` function, is the simplest way to construct
# a grid, more complicated
# (unstructured) grids can be constructed using a dictionary containing the
# vertices and the element connectivity or by for example providing a
# `gmsh` file. This will be discussed in later section, see for example
# here for [a 3d gmsh file](othergrids_nb.ipynb)
# or here [for a simple example using a dictionary](dune-corepy_nb.ipynb#An-unstructured-grid).
#
# .. index::
#    pair: Grids; Plotting
#
# %%
from matplotlib import pyplot as plt
import numpy as np

from dune.grid import structuredGrid as leafGridView
gridView = leafGridView([0, 0], [1, 1], [4, 4])
gridView.plot()

# %% [markdown]
# Often the initial grid is a bit coarse - it is very easy to refine it
# globally.
# We will discuss [local refinement](gridviews.rst#Dynamic-Local-Grid-Refinement-and-Coarsening)
# in a later section.
# Note that the ``gridView`` we have constructed can not be changed
# directly, we need to access the underlying ``hierarchicalGrid``:

# %%
from dune.fem import globalRefine
globalRefine(1,gridView.hierarchicalGrid) # the first argument is the number of refinements to perform
gridView.plot()

# %% [markdown]
# .. index::
#    pair: Functions; Grid Functions
#
# ## Grid Functions
#
# .. note:: **Grid functions** are functions that are defined
# over a given grid and are evaluated by using an element
# of the grid and local coordinate within that element.
#
# For example:
# ```
#    value = gridFunction(element,localCoordinate)
# ```
# Alternatively one can obtain a `LocalFunction` from a grid function which
# can be bound to an element and then evaluate via local coordinate:
# ```
#    localFunction = gridFunction.localFunction()
#    for e in grid.elements:
#        localFunction.bind(e)
#        value = localFunction(x)
#        localFunction.unbind()
# ```
# There are multiple ways to construct grid functions. The easiest way it
# to use UFL expression. Many methods expecting grid functions are argument
# can also directly handle UFL expression.
# We can for example integrate a UFL expression over the grid:

# %%
from ufl import SpatialCoordinate, triangle
x = SpatialCoordinate(triangle)

exact = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

from dune.fem.function import integrate
print( integrate(gridView, exact, order=5) )

# %% [markdown]
# .. index::
#    pair: I/O; VTK files

# and plot them using matplotlib or write a vtk file for postprocessing

# %%
from dune.fem.plotting import plotPointData as plot
plot(exact, grid=gridView)
gridView.writeVTK('exact', pointdata={'exact': exact})

# %% [markdown]
# In some cases it can be necessary to convert a UFL expression into a grid
# function explicitly - for example to be able to evaluate it over each
# element in the grid ([in a later section](dune-corepy_nb.ipynb)
# we provide more detail on how to access the grid interface).

# %%
from dune.fem.function import uflFunction
exact_gf = uflFunction(gridView, name="ufl", order=1, ufl=exact)
mass = 0
for element in gridView.elements:
    mass += exact_gf(element,[0.5,0.5]) * element.geometry.volume
print(mass)

# %% [markdown]
# .. index::
#    triple: Functions; Grid Functions; Decorator
#
# Another way to obtain a grid function is to use the `gridFunction`
# decorator. This can be obtained from `dune.grid` but then without UFL
# support. Using the decorator from `dune.fem.function` the resulting grid
# function can be used seamlessly within UFL expressions:

# %%
from dune.fem.function import gridFunction
@gridFunction(gridView,name="callback",order=1)
def exactLocal(element,xLocal):
    x = element.geometry.toGlobal(xLocal)
    return 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

# %% [markdown]
# we can use the same approach but with a function using global
# coordinates but can then be used like any other grid function:

# %%
@gridFunction(gridView,name="callback",order=1)
def exactGlobal(x):
    return 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

lf = exactGlobal.localFunction()

mass = 0
for element in gridView.elements:
    lf.bind(element)
    mass += lf([0.5,0.5]) * element.geometry.volume
    lf.unbind()
print(mass)

# %% [markdown]
# .. index::
#    pair: Functions; UFL functions
#
# As pointed out the `dune.fem` grid function can be used like any other
# UFL coefficient to form UFL expressions:

# %%
print( integrate(gridView, abs(exact-exactLocal), order=5) )
gf = uflFunction(gridView, name="ufl", order=1, ufl=exact+exactLocal*exactGlobal)
fig = plt.figure(figsize=(20,10))
gf.plot(figure=(fig,121))
exactGlobal.plot(figure=(fig,122))

# %% [markdown]
# Converting UFL expressions to grid functions leads to JIT code generation
# and is therefore efficient when used in other C++ algorithm (like
# `integrate`). But more complex control structures are hard to express
# within UFL. On the other hand using the `gridFunction` decorator leads
# to a callback into Python for each evaluation and is therefore much less
# efficient. An alternative approach is based on writing small C++ snippets
# implementing the grid function as described in the section on
# generating
# [C++ Based Grid Functions](cppfunctions_nb.ipynb).

# %% [markdown]
# .. index::
#    pair: Functions; Discrete Functions
#
# .. index::
#    pair: Spaces; Lagrange Space
#
# ## Discrete Functions
#
# The grid functions set up so far did not involve any
# discretization, they are exactly evaluated at the given point.
#
# .. note::
# A special type of grid function are **discrete functions** living in a
# discrete (finite dimensional) space.

# %%
from dune.fem.space import lagrange as solutionSpace
space = solutionSpace(gridView, order=2)

# %% [markdown]
# The easiest way to construct a discrete function is to use the interpolate
# method on the discrete function space.

# %%
u_h = space.interpolate(exact, name='u_h')

# %% [markdown]
# On an existing discrete function the `interpolate` method can be used to
# reinitialize it

# %%
u_h.interpolate( lambda x: 1/2.*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1 )
u_h.plot()

# %% [markdown]
# Note that in the last example the Python lambda is used as a callback
# automatically using the same concept used in the `gridFunction`
# decorator. As pointed out above there are some methods where these
# conversions are implicit and no explicit generation of a grid function
# has to be carried out.
#
# If a discrete function is already available it is possible to call `copy`
# to obtain further discrete functions:

# %%
u_h_n = u_h.copy(name="previous")

# %% [markdown]
# Finally, `clear` can be called on a discrete function which sets all
# coefficient to zero and `assign` can be used to copy all coefficients
# between two discrete function over the same space:

# %%
u_h_n.clear()
u_h_n.assign( u_h )

# %% [markdown]
# All the things we did above with grid functions can be done with discrete
# functions, e.g., evaluate locally

# %%
localUh = u_h.localFunction()
mass = 0
for element in gridView.elements:
    localUh.bind(element) # using u_h(element,[0.5,0.5]) also works
    mass += localUh([0.5,0.5]) * element.geometry.volume
    localUh.unbind()
print(mass)

# %% [markdown]
# .. index::
#    pair: Functions; Plotting
#
# .. index::
#    pair: I/O; VTK files
#
# or plot using matplotlib and write a vtk file for postprocessing
# (using binary data format to reduce size)

# %%
u_h.plot(gridLines="white")
from dune.grid import OutputType
gridView.writeVTK('uh', pointdata=[u_h])

# %% [markdown]
# Note: the discrete function `u_h` already has a `name` attribute given in
# the `interpolate` call. This is used by default in the vtk file. An
# alternative name can be given by using a dictionary as shown previously.
#
# Of course a discrete function can also be used as a coefficient in a UFL expression

# %%
print( integrate(gridView, abs(exact-u_h), order=5) )

# %% [markdown]
# The main difference between a grid function and a discrete function is
# that the latter has a `dofVector` attached to it. We can iterate over
# that easily
# %%
for d in u_h.dofVector:
    print(d,end=", ")
print()

# %% [markdown]
# In addition to `clear` and `assign` we can also add (scalar
# multiples of) discrete functions together or multiply them by a scalar
# which will directly change the underlying dofVector - this can only be
# done __in place__, so `u_h *= 2` changes the dof vector while
# `u_h*2` results in an UFL expression. Here are some examples:
# %%
u_h_n *= 2
u_h   += u_h_n
# this is equivalent to
u_h.axpy(2,u_h_n)

# %% [markdown]
# We started with `u_h_n=u_h` so after these changes each dof in
# `u_h` should be seven times its original value
# %%
for d in u_h.dofVector:
    print(d,end=", ")
print()

# %% [markdown]
# [In a later section](solvers_nb.ipynb)
# we will see how to extract the underlying dof vector in the form of a `numpy` or a `petsc` vector.

# %% [markdown]
# ## Handling constant parameters
# .. index::
#    pair: Models; Constants
#
# We have already seen how grid function can be used within UFL expressions
# and forms - in the [next section](dune-fempy_nb.ipynb) we will give another example for this in the
# context of a time dependent problems. In addition we can also use the
# `Constant` class to add constants to UFL expressions/forms which then can
# be changed without requiring any recompilation of the model. An example
# would again be in a time dependent problem a time varying coefficient.
# Being able to change the value of the time in the model without
# recompilation is crucial for an efficient code. We will demonstrate this
# here by considering
# $$ -\triangle u + \alpha u = f $$
#
# where $\alpha$ is a real valued parameter.
# We choose $f,g_N,g_D$ so that the exact solution for $\alpha=0$
# is given by
# \begin{align*}
# u(x) = \left(\frac{1}{2}(x^2 + y^2) - \frac{1}{3}(x^3 - y^3)\right) + 1
# \end{align*}
# which is the grid function ``exact`` we have already defined above.
#
# Although the problem is linear and we could just use
# ``dune.fem.assemble`` and Scipy to solve, we use a ``scheme`` instead:

# %%
from ufl import (TestFunction, TrialFunction, SpatialCoordinate, FacetNormal,
                 dx, ds, grad, div, inner, dot, sin, cos, pi, conditional)
from dune.ufl import Constant, DirichletBC
from dune.fem.scheme import galerkin as solutionScheme

alpha = Constant(1,name="mass")   # start with alpha=1
x     = SpatialCoordinate(space)
u     = TrialFunction(space)
v     = TestFunction(space)
f     = -div( grad(exact) )
g_N   = grad(exact)
n     = FacetNormal(space)
b     = f*v*dx + dot(g_N,n)*conditional(x[0]>=1e-8,1,0)*v*ds
a     = dot(grad(u), grad(v)) * dx + alpha*u*v * dx
dbc   = DirichletBC(space,exact,x[0]<=1e-8)

scheme = solutionScheme([a == b, dbc], solver='cg')
scheme.solve(target = u_h)

h1error = dot(grad(u_h - exact), grad(u_h - exact))
error = np.sqrt(integrate(gridView, h1error, order=5))
print("Number of elements:",gridView.size(0),
      "number of dofs:",space.size,"H^1 error:", error)

# %% [markdown]
#
# Note that $\alpha=1$ so we do not have the correct exact solution.
# We can print the value of a `Constant` with name `foo` via
# `scheme.model.foo` and change it's value using the same attribute:

# %%
print(scheme.model.mass)
scheme.model.mass = 0      # change to the problem which matches our exact solution
print(scheme.model.mass)

scheme.solve(target = u_h)

h1error = dot(grad(u_h - exact), grad(u_h - exact))
error = np.sqrt(integrate(gridView, h1error, order=5))
print("Number of elements:",gridView.size(0),
      "number of dofs:",space.size,"H^1 error:", error)

# %% [markdown]
# .. index:: Equations; Non-linear parabolic
#
# ## A Nonlinear Time Dependent Problem
#
# As a second example we will study the Forchheimer problem
# <cite data-cite="Kieu"></cite>
# which is a scalar, nonlinear parabolic equation
# \begin{equation}
# \partial_tu - \nabla\cdot K(\nabla u)\nabla u = f
# \end{equation}
# where the diffusion tensor is given by
# \begin{equation}
# K(\nabla u) = \frac{2}{1+\sqrt{1+4|\nabla u|}}
# \end{equation}
# and $f=f(x,t)$ is some time dependent forcing term.
# On the boundary we prescribe Neumann boundary as before
# and initial conditions $u=u_0$.
#
# We will solve this problem in variational form and using Crank Nicolson in time
#
# \begin{equation}
# \begin{split}
# \int_{\Omega} \frac{u^{n+1}-u^n}{\Delta t} \varphi
# + \frac{1}{2}K(\nabla u^{n+1}) \nabla u^{n+1} \cdot \nabla \varphi \
# + \frac{1}{2}K(\nabla u^n) \nabla u^n \cdot \nabla \varphi v\ dx \\
# - \int_{\Omega} \frac{1}{2}(f(x,t^n)+f(x,t^n+\Delta t) \varphi\ dx
# - \int_{\partial \Omega} \frac{1}{2}(g(x,t^n)+g(x,t^n+\Delta t)) v\ ds
# = 0.
# \end{split}
# \end{equation}
#
# on a domain $\Omega=[0,1]^2$. We choose $f,g$ so that the exact solution
# is given by
# \begin{align*}
# u(x,t) = e^{-2t}\left(\frac{1}{2}(x^2 + y^2) - \frac{1}{3}(x^3 - y^3)\right) + 1
# \end{align*}
# We solve this problem using a quadratic finite element space.
# We first setup this space, define the initial data,
# and a second discrete function to store
# the solution at the current time step (``u_h`` will be used to store the
# next time step):

# %%
space = solutionSpace(gridView, order=2)
u_h   = space.interpolate(0, name='u_h')
u_h_n = u_h.copy(name="previous")

x = SpatialCoordinate(space)
initial = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1

# %% [markdown]
# Next we setup the form $a=a(u,v)$ which depends on time $t$ and the
# chosen time step $\Delta t$. Since we want to be able to change these
# during the simulation (we need to do this for $t$ and might want to do
# this for $\Delta t$ we use the ``dune.ufl.Constant`` which stores a
# mutable scalar and can be used within a ufl expression:

# %%
from ufl import FacetNormal, ds, div, dot, sqrt, exp
from dune.ufl import Constant
dt = Constant(0, name="dt")    # time step
t  = Constant(0, name="t")     # current time
u = TrialFunction(space)
v = TestFunction(space)

abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2/(1 + sqrt(1 + 4*abs_du(u)))
a = ( dot((u - u_h_n)/dt, v) \
    + 0.5*dot(K(u)*grad(u), grad(v)) \
    + 0.5*dot(K(u_h_n)*grad(u_h_n), grad(v)) ) * dx

# %% [markdown]
# To define the right hand side we use the given function $u$ to define
# the forcing $f$ and the boundary flux $g$

# %%
exact = lambda t: exp(-2*t)*(initial - 1) + 1
f     = lambda s: -2*exp(-2*s)*(initial - 1) - div( K(exact(s))*grad(exact(s)) )
g = lambda s: K(exact(s))*grad(exact(s))
n = FacetNormal(space)
b = 0.5*(f(t)+f(t+dt))*v*dx + 0.5*dot(g(t)+g(t+dt),n)*v*ds

# %% [markdown]
# With the model described as a ufl form, we can again construct a scheme class
# that provides the solve method which we can use to evolve the solution from
# one time step to the next. As already mentioned the solve methods always assumes
# that the problem is non-linear and use a Newton method to solve the problem:

# %%
scheme = solutionScheme(a == b, solver='cg')

# %% [markdown]
#
# Optionally, we may want to increase the default quadrature orders which are
# `2 * space.order` for element integrals and `2 * space.order + 1` for surface
# integrals. Depending on the data this might not be enough. Then we simply
# set the integration orders by hand like in the following example, by calling
# the method setQuadratureOrders(interiorOrder, surfaceOrder).

# %%
#scheme.setQuadratureOrders( 2*space.order, 2*space.order+1 )

# %% [markdown]
#
# .. index:: Experimental Order of Convergence (EOC)
#
# Since we have forced the system towards a given solution, we can compute
# the discretization error. First we define ufl expressions for the $L^2$
# and $H^1$ norms and will use those to compute the experimental order of
# convergence of the scheme by computing the time evolution on different grid
# levels. We use ``dune.fem.function.uflFunction`` to generate code for the
# errors which can be more efficient than using the ufl expressions
# directly.

# %%
endTime    = 0.25
exact_end  = exact(endTime)
l2error = uflFunction(gridView, name="l2error", order=u_h.space.order,
                      ufl=dot(u_h - exact_end, u_h - exact_end))
h1error = uflFunction(gridView, name="h1error", order=u_h.space.order,
                      ufl=dot(grad(u_h - exact_end), grad(u_h - exact_end)))

# %% [markdown]
# .. index::
#    pair: I/O; VTK files
#
# Now we evolve the solution from time $t=0$ to $t=T$ in a loop. Since the
# PDE model has time dependent coefficient (through the forcing term), we
# need to update the `t` constant used to define the model before each step.
# A second constant we used to define the model was `dt` which defines the
# time step. We keep this constant and set it to $0.002$ at the beginning
# of the simulation. This value could be changed in each time step:

# %%
scheme.model.dt = 0.01
time = 0
u_h.interpolate(initial)
while time < (endTime - 1e-6):
    scheme.model.t = time
    u_h_n.assign(u_h)
    scheme.solve(target=u_h)
    time += scheme.model.dt

errors = [np.sqrt(e) for e in integrate(gridView, [l2error,h1error], order=5)]
print('grid size:', gridView.size(0))
print('\t | u_h - u | =', '{:0.5e}'.format(errors[0]))
print('\t | grad(uh - u) | =', '{:0.5e}'.format(errors[1]))
u_h.plot()
gridView.writeVTK('forchheimer', pointdata={'u': u_h, 'l2error': l2error, 'h1error': h1error})

# %% [markdown]
# We can refine the grid once and recompute the solution on the finer grid.

# %%
globalRefine(1, gridView.hierarchicalGrid)
scheme.model.dt /= 2
u_h.interpolate(initial)
time = 0
while time < (endTime - 1e-6):
    scheme.model.t = time
    u_h_n.assign(u_h)
    scheme.solve(target=u_h)
    time += scheme.model.dt

errorsFine = [np.sqrt(e) for e in integrate(gridView, [l2error,h1error], order=5)]
print('grid size:', gridView.size(0))
print('\t | u_h - u | =', '{:0.5e}'.format(errorsFine[0]))
print('\t | grad(uh - u) | =', '{:0.5e}'.format(errorsFine[1]))
u_h.plot()

# %% [markdown]
# To check that everything is working as expected let's compute the
# **experimental order of convergence** (convergence should be cubic for
# $L^2$-error and quadratic for the $H^1$-error:

# %%
eocs = [ round(np.log(fine/coarse)/np.log(0.5),2)
         for fine,coarse in zip(errorsFine,errors) ]
print("EOCs:",eocs)

# %% [markdown]
# .. index:: Grids; Available Grids
#
# ## Listing Available Dune Components
# The available realization of a given interface, i.e., the available
# grid implementations, depends on the modules found during configuration.
# Getting access to all available components is straightforward:

# %%
from dune.utility import components
# to get a list of all available components:
components()
# to get for example all available grid implementations:
components("grid")
# %%

# %% [markdown]
# .. index:: Spaces; Available Spaces
#
# Available discrete function spaces are:

# %%
components("space")

# %% [markdown]
# .. index:: Functions; Available Functions
#
# Available discrete functions are:

# %%
components("function")

# %% [markdown]
# .. index:: Schemes; Available Schemes
#
# Available schemes:

# %%
components("scheme")
