# %% [markdown]
#
# # Tweaking the (non-) linear solvers
# Here we look at how to tweak the internal linear and nonlinear solvers.
# The `solve` method on the `schemes` uses a Newton method to solve
# the (non)linear PDE model. Parameters set during the construction of a `scheme`
# can be used to customize the Newton solver and the underlying solver for
# the linear problem in each step, i.e., setting line-search parameters or
# preconditioners. A number of linear solvers and some simple
# preconditioners are directly available and others can be used through
# additional software packages. This is done by changing the
# linear algebra backend (e.g. using `petsc`) by changing using the
# `storage` parameter during construction of the discrete space.
# The simplicity of the default backend (`numpy`) allows the usage
# Numpy/Scipy data structures on the python side. A degrees of freedom
# vector (dof vector) can be retrieved from a discrete function with the
# default storage by using the `as_numpy` method. Similar methods are available
# for the other storages, i.e., `as_istl,as_petsc`. The same methods are
# also available to retrieve the underlying matrix structures of linear
# operators. Here we will discuss the `as_numpy` attribute and how to use
# this to implement part of the algorithm for solving PDEs on the Python side.
# A short description on the other storages follows at the end of this
# section, with more detail [to follow](solvers-external_nb.ipynb).
#
#
# We will mostly revisit the nonlinear time dependent problem studied at the end of the
# [concepts section](concepts_nb.ipynb) which after discretizing in time had the variational formulation
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
# on a domain $\Omega=[0,1]^2$. We choose $f,g$ so that the exact solution
# is given by
# \begin{align*}
# u(x,t) = e^{-2t}\left(\frac{1}{2}(x^2 + y^2) - \frac{1}{3}(x^3 - y^3)\right) + 1
# \end{align*}
# The following code was described in the [concepts section](concepts_nb.ipynb)

# %%
import numpy, sys, io
import matplotlib.pyplot as plt

from dune.grid import structuredGrid as leafGridView
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import gridFunction
from dune.fem import integrate
from dune.ufl import Constant
from ufl import TestFunction, TrialFunction, SpatialCoordinate, FacetNormal, \
                dx, ds, div, grad, dot, inner, sqrt, exp, sin,\
                conditional

gridView = leafGridView([0, 0], [1, 1], [12, 12])
space = lagrange(gridView, order=2)

x = SpatialCoordinate(space)
initial = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
exact   = lambda t: exp(-2*t)*(initial - 1) + 1

u_h   = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")

u = TrialFunction(space)
v = TestFunction(space)
dt = Constant(1, name="dt")       # time step
t  = Constant(0, name="time")     # current time

abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2/(1 + sqrt(1 + 4*abs_du(u)))
a = ( dot((u - u_h_n)/dt, v) \
    + 0.5*dot(K(u)*grad(u), grad(v)) \
    + 0.5*dot(K(u_h_n)*grad(u_h_n), grad(v)) ) * dx

f = lambda s: -2*exp(-2*s)*(initial - 1) - div( K(exact(s))*grad(exact(s)) )
g = lambda s: K(exact(s))*grad(exact(s))
bf = 0.5*(f(t)+f(t+dt))*v*dx

n = FacetNormal(space)
bg = 0.5*dot(g(t)+g(t+dt),n)*v*ds

# %% [markdown]
# When creating a scheme, it is possible to set the linear solver as well as
# parameters for the internal Newton solver and the linear solver
# and preconditioning. See a list of available solvers and preconditioning
# methods at the end of this section. The default solver is a `gmres`
# method but since this problem is symmetric a `cg` solver is more
# efficient. In addition we will use a Jacobi preconditioner. We also set
# the tolerance for the Newton method to $10^{-9}$:

# %%
scheme = galerkin(a == bf+bg, solver='cg',
                  parameters={"newton.linear.preconditioner":"jacobi",
                              "newton.tolerance":1e-8},
                 )

endTime    = 0.25
exact_end  = exact(endTime)
l2error = gridFunction(name="l2error", expr=dot(u_h - exact_end, u_h - exact_end))
h1error = gridFunction(name="h1error", expr=dot(grad(u_h - exact_end), grad(u_h - exact_end)))

# %% [markdown]
# We define a function to evolve the solution from time 0 to the end time.
# The first argument is a class with a `solve` method that moves the
# solution from one time level to the next - i.e., solves the non-linear
# problem for $u^{n+1}$ given $u^n$:

# %%
import time as timer
def evolve(scheme, u_h, u_h_n, endTime):
    start = timer.time()
    time = 0
    nIter,linIter = 0,0
    steps = 0
    while time < (endTime - 1e-6):
        scheme.model.time = time
        t.value = time
        # assert t.value == time
        u_h_n.assign(u_h)
        info = scheme.solve(target=u_h)
        nIter += info["iterations"]
        linIter += info["linear_iterations"]
        time += scheme.model.dt
        steps += 1
    print(f"steps={[steps,scheme.model.dt]}, runtime={timer.time()- start}, iterations={[nIter,linIter]}",flush=True)

# %% [markdown]
# We can simply use the `galerkin' scheme instance in this function to produce the solution
# at the final time this is comparable to the code shown previously in
# [concepts section](concepts_nb.ipynb).
# We combine this with a loop to compute the error over
# two grids and estimate the convergence rate:

# %%
dt.value = 0.05

errors = 0,0
loops = 2
for eocLoop in range(loops):
    t.value = 0
    u_h.interpolate(initial)
    evolve(scheme, u_h, u_h_n, endTime)
    errors_old = errors
    errors = [sqrt(e) for e in integrate([l2error,h1error])]
    if eocLoop == 0:
        eocs = ['-','-']
    else:
        eocs = [ round(numpy.log(e/e_old)/numpy.log(0.5),2) \
                 for e,e_old in zip(errors,errors_old) ]
    print('Forchheimer: step:', eocLoop, ', size:', gridView.size(0),flush=True)
    print('\t | u_h - u | =', '{:0.5e}'.format(errors[0]), ', eoc =', eocs[0], flush=True)
    print('\t | grad(uh - u) | =', '{:0.5e}'.format(errors[1]), ', eoc =', eocs[1])
    u_h.plot()
    if eocLoop < loops-1:
        gridView.hierarchicalGrid.globalRefine(1)
        dt.value /= 2

# %% [markdown]
#
# ## Python sided Newton solver
# We can also only use the internal linear solvers and write the Newton
# method in Python - here is a very simple version of that.
# The central new method is `linear` available on the scheme that returns a
# linearized version of the model around a given grid function `ubar`.
# Much more details on the methods used are given in the next section.

# %%
class Scheme:
    def __init__(self, scheme, u0):
        self.scheme = scheme
        self.model = scheme.model
        self.jacobian = scheme.linear(ubar=u0) # obtain a linear operator for the Newton method

    def solve(self, target):
        # create a copy of target for the residual
        res = target.copy(name="residual")
        dh  = target.copy(name="direction")

        n = 0
        linIter = 0
        while True:
            self.scheme(target, res)
            absF = res.scalarProductDofs(res)
            if absF < 1e-8**2: # this is the same tolerance we set above for the built-in Newton solver
                break
            self.scheme.jacobian(target,self.jacobian)  # assemble the linearization
            dh.clear()                                  # zero is good initial guess for iterator solver
            info = self.jacobian.solve(target=dh,rhs=res)
            linIter += info["linear_iterations"]
            target -= dh
            n += 1
        return {"iterations":n, "linear_iterations":linIter}

t.value = 0
u_h.interpolate(initial)                # reset u_h to initial
scheme_cls = Scheme(scheme,u0=u_h)

evolve(scheme_cls, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(Newton 1) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
  *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]), flush=True)

# %% [markdown]
# A second version using a linearization of the non-linear problem
# including an affine shift, i.e., the first two terms in the Taylor
# expansion around a given $\bar{u}$:
# $$ L[u] = DS[\bar{u}](u-\bar{u}) + S[\bar{u}] $$

# %%
from dune.fem.scheme import linearized
class Scheme:
    def __init__(self, scheme, u0):
        self.scheme = scheme
        self.model = scheme.model
        # let's use a `gmres` solver for demonstration purposes
        self.linear = linearized(scheme, ubar=u0,
                      parameters={ "method":"gmres", "gmres.restart":50 })
        self.linear.setup(u0)

    def solve(self, target):
        prev = target.copy()
        n = 0
        linIter = 0
        # we are not computing the residual directly so need a different
        # stopping criteria, i.e., size of the step dh
        while True:
            # setup L[u] = DS[u^n](u-u^n) + S[u^n]
            # Then L[u] = 0 leads to
            # u = (DS[u^n])^{-1} (DS[u^n]u^n-S[u^n])
            #   = u^n - DS[u^n])^{-1}S[u^n]
            #   = u^{n+1}
            self.linear.setup(target)
            info = self.linear.solve(target=target)
            linIter += info["linear_iterations"]
            n += 1

            # subtract the new iterate from the previous one
            prev.sub( target )
            absDiff = prev.scalarProductDofs(prev)
            if absDiff < 1e-8**2: # this is the default for the built-in Newton solver
                break
            prev.assign(target)
        return {"iterations":n, "linear_iterations":linIter}

t.value = 0
u_h.interpolate(initial)                # reset u_h to initial
u_h_n.assign(u_h)
scheme_cls = Scheme(scheme,u_h)

evolve(scheme_cls, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(Newton 2) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
  *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]),flush=True)

# %% [markdown]
# __Note__: the parameters used to tweak the linear solver set when constructing
# the `scheme` are passed on to the linear solver when calling
# `scheme.linear` but can be overwritten in the call of `scheme.linear`.
# So to obtain a linearization (around zero by default) using a `gmres`
# method when calling `jac.solve` we could call

# %%
model = ( inner(grad(u), grad(v)) + (u - conditional(dot(x,x)<0.6,10.,0.))*v ) * dx
# verbose output for the linear solver is too long so use default `False` setting
sc  = galerkin(model==0, solver="cg",
               parameters={"newton.linear.verbose":False,"newton.verbose":True}
              )
sol = space.interpolate(0, name="u_h")
jac = linearized(sc, parameters={ "method":"gmres", "gmres.restart":100, })
info = jac.solve(target=sol)
print(info,flush=True)

# %% [markdown]
# Since we setup the `scheme` using `solver="cg"` calling the `solve`
# method will use the Newton method + a `cg` solver.
#
# __Note__: since the problem is linear the Newton method produces no real
# overhead.

# %%
sol.clear()
info = sc.solve(target=sol)
print(info,flush=True)

# %% [markdown]
# .. index::
#    triple: Solvers; Dirichlet; Conditions
#
# .. index::
#    pair: Boundary; Linear Solvers
#
# ## Handling Dirichlet boundary conditions
# We first look at a simple Poisson problem with Dirichlet BCs
# to show how to use the internal solvers in this setting.
# We solve $-\triangle u=10\chi_\omega$ where $\chi_\omega$ is a characteristic
# function with $\omega=\{x\colon |x|^2<0.6\}$. For the boundary
# we prescribe trivial Neumann at the top and bottom boundaries
# and Dirichlet values $u=-1$ and $u=1$ at the left and right
# boundaries, respectively.

# %%
from dune.ufl import DirichletBC
from scipy.sparse.linalg import cg as solver
model = ( inner(grad(u), grad(v)) -
          conditional(dot(x,x)<0.6,10.,0.) * v ) * dx
dbcs  = [ DirichletBC(space,-1,1), DirichletBC(space, 1,2) ]
op  = galerkin([model==0, *dbcs], space)
sol = space.interpolate(0, name="u_h")
lin = op.linear()

# %% [markdown]
# So far everything is as before. Dirichlet boundary conditions
# are handled in the matrix through changing all rows
# associated with boundary degrees of freedom to unit rows -
# associated columns are not changed so the matrix will not be symmetric anymore.
# For solving the system we need to modify the right hand side
# and the initial guess for the iterative solver to include
# the boundary values (to counter the missing symmetry).
# We can use the first of the three versions of the
# `setConstraints` methods on the scheme class discussed in the section on
# [more general boundary conditions](boundary_nb.ipynb#Accessing-the-Dirichlet-degrees-of-freedom).

# %%
op.setConstraints(sol)

# setup the right hand side
rhs = sol.copy()
op(sol,rhs)
rhs *= -1
op.setConstraints(rhs)
lin.solve(rhs=rhs,target=sol)
sol.plot()

# %% [markdown]
# Let's revisit the first implementation of the Newton method discussed
# above but change the time dependent problem to require Dirichlet
# conditions on the whole boundary:

# %%
scheme = galerkin([a == bf, DirichletBC(space,exact(t+0.5*dt))], solver='cg')

class Scheme:
    def __init__(self, scheme, u0):
        self.scheme = scheme
        self.model = scheme.model
        self.jacobian = scheme.linear(ubar=u0)  # obtain a linear operator for the Newton method

    def solve(self, target):
        # create a copy of target for the residual
        res = target.copy(name="residual")
        dh  = target.copy(name="direction")

        n = 0
        linIter = 0
        while True:
            self.scheme.setConstraints(target)
            self.scheme(target, res)
            absF = res.scalarProductDofs(res)
            if absF < 1e-8**2: # this is the default for the built-in Newton solver
                break
            self.scheme.jacobian(target,self.jacobian)  # assemble the linearization
            dh.clear()                                  # zero is good initial guess for iterator solver
            info = self.jacobian.solve(target=dh,rhs=res)
            linIter += info["linear_iterations"]
            target -= dh
            n += 1
        return {"iterations":n, "linear_iterations":linIter}
t.value = 0
u_h.interpolate(initial)                # reset u_h to initial
scheme_cls = Scheme(scheme,u0=u_h)
evolve(scheme_cls, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(Dirichlet, Newton 1) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
  *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]),flush=True)

# %% [markdown]
# ## Other internal solvers (petsc)
# .. index::
#    pair: Solvers; Petsc
#
# The following requires that a PETSc installation was found during the
# configuration of ``dune``.

# %%
from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
try:
    assertCMakeHave("HAVE_PETSC")
    petsc = True
except ConfigurationError:
    print("Dune not configured with petsc - skipping example")
    petsc = False

# %% [markdown]
# We provide internal bindings to parts of the PETSc package and show this
# first. How to make use of `petsc4py` is shown
# [in a later section](solvers-external_nb.ipynb).
#
# Switching to a storage based on the PETSc solver package and solving the
# system using the dune-fem bindings can be achieved by using the
# ``storage`` argument to the space constructor.

# %%
if petsc:
    spacePetsc = lagrange(gridView, order=2, storage='petsc')
    # We use an sor preconditioner:
    schemePetsc = galerkin(a == bf+bg, space=spacePetsc,
                           parameters={"linear.preconditioning.method":"sor"})
    u_h = spacePetsc.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")
    evolve(schemePetsc, u_h, u_h_n, endTime)
    error = u_h - exact_end
    print("Forchheimer(petsc) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
      *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))

# %% [markdown]
# ## Other internal solvers (istl)
# .. index::
#    pair: Solvers; dune-istl
#
# We can just as easily use the `dune-istl` package to solve the system
# (`dune-istl` should be available when installing `dune-fem`):

# %%
spaceIstl = lagrange(gridView, order=2, storage='istl')
# We use an amg preconditioner:
schemeIstl = galerkin(a == bf+bg, space=spaceIstl,
                      parameters={"linear.preconditioning.method":"amg-jacobi"})
u_h = spaceIstl.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")
evolve(schemeIstl, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(istl) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
      *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))




# %% [markdown]
# .. index::
#    triple: I/O; Logging; Parameters

# ## Parameters for fine-tuning the internal solvers
# Sometimes it is necessary to extract which parameters were read and which
# values were used, e.g., for debugging purposes like finding spelling
# in the parameters provided to a scheme.
# Note that this information can only be reliably obtained after usage of
# the scheme, e.g., after calling solve as shown in the example below.
# To add logging to a set of parameters passed to a `scheme` one simply
# needs to add a `logging` key to the parameter dictionary provided to the scheme
# with a tag (string) that is used in the output.
#
# As an example we will solve the simple Laplace equation from the
# introduction but pass some preconditioning parameters to the scheme.

# %%
import dune.fem
from dune.grid import structuredGrid
from dune.fem.space import lagrange
from ufl import (TestFunction, TrialFunction, SpatialCoordinate,
                 dx, grad, inner, dot, sin, cos, pi )
gridView = structuredGrid([0, 0], [1, 1], [200, 200])
space = lagrange(gridView, order=1)
u_h   = space.interpolate(0, name='u_h')
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)

f = (8*pi**2+1) * cos(2*pi*x[0])*cos(2*pi*x[1])
a = ( inner(grad(u),grad(v)) + u*v ) * dx
l = f*v * dx
scheme = galerkin( a==l, solver="cg", parameters=
                   {"newton.linear.tolerance": 1e-12,
                    "newton.linear.verbose": True,
                    "newton.linear.preconditioning.method": "jacobi",
                    "newton.linear.errormeasure": "relative",
                    "logging": "log-jacobi"
                   } )
info = scheme.solve(target=u_h)

# %% [markdown]
# We use the `pprint` (pretty print) module if available to get nicer
# output.

# %%
try:
    from pprint import pprint as _pprint
    pprint = lambda *args,**kwargs: _pprint(*args,**kwargs,width=200,compact=False)
except ImportError:
    pprint = print

pprint(dune.fem.parameter.log())

# %% [markdown]
# Note above that all parameters are printed including some default ones
# used in other parts of the code. If multiple schemes with different
# `logging` parameter strings are used, all would be shown using the `log`
# method as shown above.
# To access only the parameters used in the scheme simply use
# either `dune.fem.parameter.log()["tag"])` or access the parameter log
# through the scheme:

# %%
pprint(scheme.parameterLog())

# %% [markdown]
# One can easily reuse these parameters to construct another scheme by
# converting the result of the above call to a dictionary.
# As an example change the above problem to a PDE with Dirichlet conditions
# but turn off verbose output of the solver.
#
# .. note:: the `logging` parameter has to be set if we want to use the
# `parameterLog` method on the scheme.

# %%
param = dict(scheme.parameterLog()) # this method returns a set of pairs which we can convert to a dictionary
param["logging"] = "Dirichlet" # only needed to use the `parameterLog` method
param["newton.linear.verbose"] = False
scheme2 = galerkin( [a==l,DirichletBC(space,0)], parameters=param )
u_h.clear()
info = scheme2.solve(target=u_h)
pprint(scheme2.parameterLog())

# %% [markdown]
# .. tip:: To get information about available values for some parameters
# (those with string arguments) a possible approach is to provide a non valid
# string, e.g., `"help"`.

# %%
try:
    scheme = galerkin( a==l, solver="cg", parameters=
                       {"newton.linear.tolerance": 1e-12,
                        "newton.linear.verbose": True,
                        "newton.linear.preconditioning.method": "help",
                        "newton.linear.errormeasure": "relative",
                        "logging": "precon"
                       } )
    scheme.solve(target=u_h)
except RuntimeError as rte:
    print(rte)

# %% [markdown]
# .. index::
#    triple: Solvers; Available solvers; Parameters
#
# ## Available solvers and parameters
# Upon creation of a discrete function space one also have to specifies the
# storage which is tied to the solver backend.
# As mentioned, different linear algebra backends can be accessed by changing setting the
# `storage` parameter during construction of the discrete space. All
# discrete functions and operators/schemes based on this space will then
# use this backend. Available backends are `numpy,istl,petsc`.
# Note that not all methods which are available in `dune-istl` or `PETSc` have been forwarded
# to be used with `dune-fem`.

# %%
space = lagrange(gridView, order=2)

# %% [markdown]
# Switching is as simple as passing `storage='istl'` or  `storage='petsc'`.
# Here is a summary of the available backends
#
# | Solver    | Description |
# | --------- | ------------------------------------------------------------------- |
# | numpy     | the storage is based on a raw C pointer which can be                |
# |           | directly accessed as a numpy.array using the Python buffer protocol |
# |           | To change the underlying vector of a discrete function 'u_h' use    |
# |           | 'uh.as_numpy[:]'.                                                   |
# |           | As shown in the examples, linear operators return a                 |
# |           | `scipy.sparse.csr_matrix` through the 'as_numpy' property.          |
# | istl      | data is stored in a block vector/matrix from the dune.istl package  |
# |           | Access through 'as_istl'                                            |
# | petsc     | data is stored in a petsc vector/matrix which can also be used with |
# |           | petsc4py on the python side using 'as_petsc'                        |

# %% [markdown]
# .. index:: Solvers; Parameters

# When creating a scheme, there is the possibility to select a linear
# solver for the internal Newton method.
# In addition the behavior of the solver can be customized through a
# parameter dictionary. This allows to set tolerances, verbosity, but also
# which preconditioner to use.
#
# For details see the help available for a scheme:

# %%
help(scheme)
