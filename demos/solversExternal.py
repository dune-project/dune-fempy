# %% [markdown]
#
# # Alternative Solvers (Scipy and Petsc)
# Here we look at different ways of solving PDEs using external
# packages.
# Different linear algebra backends can be accessed by changing setting the
# `storage` parameter during construction of the discrete space. All
# discrete functions and operators/schemes based on this space will then
# use this backend. Available backends are `numpy,istl,petsc`. The default is
# `numpy` which uses simple data structures and linear solvers implemented in
# the `dune-fem` package - these backends were discussed
# [in the previously section](solversInternal_nb.ipynb).
#
# As discussed there, a degrees of freedom
# vector (dof vector) can be retrieved from a discrete function over the
# `numpy` space by using the `as_numpy` method. Similar methods are available
# for the other storages, i.e., `as_istl,as_petsc`. The same methods are
# also available to retrieve the underlying matrix structures of linear
# operators.
#
# We start with the same setup as used [in the previously section](solversInternal_nb.ipynb)
# (provided in script [solversIntro.py](solversIntro.py)):

# %%
from solversIntro import *

# %% [markdown]
# Having defined the model, we can now construct an `operator`. We will
# use the `galerkin` operator from `dune.fem.operator` instead of the
# `scheme` used so far, since we do not need the `solve` method since we
# want to use external solvers in this section:

# %%
from dune.fem.operator import galerkin
operator = galerkin(a - (bf+bg))

# %% [markdown]
# ## Accessing underlying data structures
#
# The most important step is accessing the data structures setup on the C++
# side in Python. In this case we would like to use the underlying dof vector from
# a discrete function as numpy arrays and system matrices assembled by the
# schemes and operators as scipy sparse matrices.
# In the [introduction](dune-fempy_nb.ipynb) we already discussed the `as_numpy` method.
# So

# %%
vecu_h = u_h.as_numpy

# %% [markdown]
# provides access to the underlying dof vector without copying. So changes
# to the numpy array `vecu_h` carries over to the discrete function `u_h`.
# Just remember to make changes using `vecu_h[:]` to change the actual
# memory buffer.
#
# An instance of `dune.fem.operator` describing an operator `L`
# provides a method `linear` which returns an object that
# stores a sparse matrix structure. The object describes the operator
# linearized around zero.
# To linearize around a different value use the `jacobian` method on the
# `operator` that linearized
# the the operator `L` around a given grid function `ubar` and fills the
# linear operator structure passed in as second argument. It is also
# possible to pass ``assemble=False`` to the ``linear`` method to avoid an
# the linearization around zero to reduce computational cost:

# %%
linOp = operator.linear()                  # linearized around 0
linOp = operator.linear(assemble=False)    # empty (non valid) linear operator
operator.jacobian(space.zero, linOp)       # linearized around zero

# %% [markdown]
# Here we linearize around zero. But that argument could be any grid
# function. A second version of this method will return an addition
# discrete function `rhs` which equals `-L[ubar]` such that `DL[ubar](u-ubar) - rhs` are the first
# terms in the Taylor expansion of the operator `L`:

# %%
rhs = u_h.copy()
operator.jacobian(u_h, linOp, rhs)

# %% [markdown]
# One can now easily access the underlying sparse matrix by again using
# `as_numpy` (and again the underlying data buffers are not copied):

# %%
A = linOp.as_numpy
print(type(A))
# plt.spy(A, precision=1e-8, markersize=1)

# %% [markdown]
# ## Defining discrete function with given external dof vector
#
# As a last step we sometimes need to construct a discrete function based
# on a given dof vector. Here a simple example where we use a `numpy`
# random vector as dof vector for a discrete function:

# %%
v_coeff = numpy.random.rand(space.size)
vh = space.function("random", dofVector=v_coeff)
vh.plot()

# %% [markdown]
# Again the dof vector is not copied so changing the `numpy` array will
# directly change the discrete function

# %%
v_coeff[range(space.size)] = range(space.size)
vh.plot()

# %% [markdown]
# Now we have all the ingredients to write a simple Newton solver to solve
# our non-linear time dependent PDE.

# %% [markdown]
# .. index::
#    pair: Solvers; Scipy
#
# ## Using Scipy
# In our first example, we implement a simple Newton Krylov solver using a linear solver from
# Scipy. We can use the `as_numpy` method to access the degrees of freedom as
# Numpy vector based on the `python buffer protocol`. So no data is copied
# and changes to the dofs made on the Python side are automatically carried
# over to the C++ side.

# %%
import numpy as np
from scipy.sparse.linalg import spsolve as solver
class Scheme:
    def __init__(self, op):
        self.model = op.model
        self.jacobian = op.linear()
        self.op = op

    def solve(self, target):
        # create a copy of target for the residual
        res = target.copy(name="residual")

        # extract numpy vectors from target and res
        sol_coeff = target.as_numpy
        res_coeff = res.as_numpy

        n = 0
        while True:
            self.op(target, res)
            absF = numpy.sqrt( np.dot(res_coeff,res_coeff) )
            if absF < 1e-6:
                break
            self.op.jacobian(target,self.jacobian)
            sol_coeff -= solver(self.jacobian.as_numpy, res_coeff)
            n += 1
        return {"iterations":n, "linear_iterations":n} # direct linear solver requires no additional iterations

scheme_cls = Scheme(operator)

u_h.interpolate(initial)                # reset u_h to initial
evolve(scheme_cls, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(numpy) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
  *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))

# %% [markdown]
# We can also use a non linear solver from the Scipy package

# %%
from scipy.optimize import newton_krylov
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import cg as solver

class Scheme2:
    def __init__(self, op):
        self.op = operator
        self.model = op.model
        self.res = u_h.copy(name="residual")
        self.n = 0

    # non linear function
    def f(self, x_coeff):
        # the following converts a given numpy array
        # into a discrete function over the given space
        x = space.function("tmp", dofVector=x_coeff)
        self.op(x, self.res)
        self.n += 1
        return self.res.as_numpy

    # class for the derivative DS of S
    class Df(LinearOperator):
        def __init__(self, op, x_coeff):
            self.op = op
            self.shape = (x_coeff.shape[0], x_coeff.shape[0])
            self.dtype = x_coeff.dtype
            x = space.function("tmp", dofVector=x_coeff)
            self.jacobian = self.op.linear()
            self.update(x_coeff,None)
        # reassemble the matrix DF(u) given a DoF vector for u
        def update(self, x_coeff, f):
            x = space.function("tmp", dofVector=x_coeff)
            self.op.jacobian(x, self.jacobian)
        # compute DS(u)^{-1}x for a given DoF vector x
        def _matvec(self, x_coeff):
            return solver(self.jacobian.as_numpy, x_coeff, rtol=1e-10, atol=1e-10)[0]

    def solve(self, target):
        sol_coeff = target.as_numpy
        # call the newton krylov solver from scipy
        sol_coeff[:] = newton_krylov(self.f, sol_coeff,
                    f_tol=1e-6,
                    inner_M=self.Df(self.op,sol_coeff))
        return {"iterations":self.n, "linear_iterations":-1} # not keeping track of linear iterations

scheme2_cls = Scheme2(operator)
u_h.interpolate(initial)
evolve(scheme2_cls, u_h, u_h_n, endTime)
error = u_h - exact_end
print("Forchheimer(scipy) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
  *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))

# %% [markdown]
# .. index::
#    pair: Solvers; Petsc4Py
#
# ## Using Petsc4Py
# The following requires that a PETSc installation was found during the
# configuration of ``dune``. Furthermore the Python package ``petsc4py`
# is required:

# %%
from dune.common.checkconfiguration import assertCMakeHave, ConfigurationError
try:
    assertCMakeHave("HAVE_PETSC")
    petsc = True
except ConfigurationError:
    print("Dune not configured with petsc - skipping example")
    petsc = False
try:
    import petsc4py
    petsc4py.init(sys.argv)
    from petsc4py import PETSc
except ModuleNotFoundError:
    print("petsc4py module not found -- skipping example")
    petsc4py = None

# %% [markdown]
# To use `petsc4py` we need to which to a storage based on the PETSc package
# using the ``storage`` argument to the space constructor.

# %%
if petsc:
    spacePetsc = lagrange(gridView, order=2, storage='petsc')
    u_h = spacePetsc.interpolate(initial, name='u_h')
    u_h_n = u_h.copy(name="previous")

# %% [markdown]
# Implementing a Newton Krylov solver using the binding provided by
# `petsc4py` is now straightforward:

# %%
if petsc4py is not None and petsc:
    class Scheme3:
        def __init__(self, scheme):
            self.model = scheme.model
            self.jacobian = scheme.linear()
            self.ksp = PETSc.KSP()
            self.ksp.create(PETSc.COMM_WORLD)
            # use conjugate gradients method
            self.ksp.setType("cg")
            # and incomplete Cholesky
            self.ksp.getPC().setType("icc")
            self.ksp.setOperators(self.jacobian.as_petsc)
            self.ksp.setFromOptions()
        def solve(self, target):
            res = target.copy(name="residual")
            sol_coeff = target.as_petsc
            res_coeff = res.as_petsc
            n = 0
            while True:
                schemePetsc(target, res)
                absF = numpy.sqrt( res_coeff.dot(res_coeff) )
                if absF < 1e-6:
                    break
                schemePetsc.jacobian(target, self.jacobian)
                self.ksp.solve(res_coeff, res_coeff)
                sol_coeff -= res_coeff
                n += 1
            return {"iterations":n, "linear_iterations":-1} # not keeping track of linear iterations

    schemePetsc = galerkin(a==bf+bg, domainSpace=spacePetsc,rangeSpace=spacePetsc)
    u_h.interpolate(initial)
    scheme3_cls = Scheme3(schemePetsc)
    evolve(scheme3_cls, u_h, u_h_n, endTime)
    error = u_h - exact_end
    print("Forchheimer(petsc4py 1) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
      *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))

# %% [markdown]
# Using the petsc4py bindings for the non linear KSP solvers from PETSc
# %%
if petsc4py is not None and petsc is not None:
    class Scheme4:
        def __init__(self, scheme):
            self.scheme = scheme
            self.model = scheme.model
            self.res = scheme.rangeSpace.interpolate([0],name="residual")
            self.jacobian = scheme.linear()
            self.snes = PETSc.SNES().create()
            self.snes.setFunction(self.f, self.res.as_petsc.duplicate())
            self.snes.setUseMF(False)
            self.snes.setJacobian(self.Df, self.jacobian.as_petsc, self.jacobian.as_petsc)
            self.snes.getKSP().setType("cg")
            self.snes.setTolerances(atol=1e-6)
            self.snes.setFromOptions()

        def f(self, snes, x, f):
            # setup discrete function using the provide petsc vectors
            inDF = self.scheme.domainSpace.function("tmp",dofVector=x)
            outDF = self.scheme.rangeSpace.function("tmp",dofVector=f)
            self.scheme(inDF,outDF)

        def Df(self, snes, x, m, b):
            inDF = self.scheme.domainSpace.function("tmp",dofVector=x)
            self.scheme.jacobian(inDF, self.jacobian)
            return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

        def solve(self, target):
            sol_coeff = target.as_petsc
            self.res.clear()
            self.snes.solve(self.res.as_petsc, sol_coeff)
            return {"iterations":self.snes.getIterationNumber(),
                    "linear_iterations":self.snes.getLinearSolveIterations()}

    u_h.interpolate(initial)
    scheme4_cls = Scheme4(schemePetsc)
    evolve(scheme4_cls, u_h, u_h_n, endTime)
    error = u_h - exact_end
    print("Forchheimer(petsc4py 2) size: ", gridView.size(0), "L^2, H^1 error:",'{:0.5e}, {:0.5e}'.format(
      *[ sqrt(e) for e in integrate([error**2,inner(grad(error),grad(error))]) ]))
