.. index:: Operator; Operator API

API for schemes and operators
=============================

Abstract concepts
-----------------

In the following we give a short mathematical description of the operator concepts used
in this module. This can be skipped and the description of the methods
further down should still be understandable.

The operator classes represent operators :math:`L\colon V_h\to W_h^*` where
the domain and range space :math:`V_h,W_h` are discrete function spaces and
:math:`W_h^*` is the dual of :math:`W_h`.
We will in general not differentiate between duals and the actual function
because the range space is finite dimensional.
In the most common example these operators are given by forms :math:`a(v_h,w_h)`
which are linear in the second argument (test function) but can be nonlinear
in the first argument (trial function).
In this case :math:`L[v_h] := a(v_h,\cdot)\in W_h^*` or for any
:math:`w_h\in W_h\colon <L[v_h],w_h> := a(v_h,w_h)`.
This is provided by ``dune.fem.operator.galerkin`` for forms defined in ufl.

Mainly operators can be evaluated, i.e., we can compute :math:`w_h:=L[v_h]` using
the ``__call__`` method: ``operator(v_h,w_h)``. In addition it is possible to
compute the linearization :math:`DL[\bar{v}_h]` around some function
:math:`\bar{v}_h\in V_h`. This is done using the ``jacobian`` method on the operator.
The linearization is now a bilinear form :math:`\bar{a}(v_h,w_h;\bar{v_h})`.

.. note:: both the ``__call__`` and the ``jacobian`` method can be called using
          grid functions instead of discrete functions :math:`v_h` and :math:`\bar{v}_h`,
          respectively.

Operator without DBC:
---------------------

In the following we focus on the operators returned by
``dune.fem.operator.galerkin``.
The minimal requirements for constructing such an operator is to provide a ufl form
where the domain and range space is extracted from the trial and test
function. In this case these have to be constructed based on discrete function
spaces $V_h,W_h$, respectively, so one can not use the abstract
``dune.ufl.Space``. It is instead possible to provide discrete spaces as extra arguments
``domainSpace,rangeSpace`` to override the spaces uses to construct the
trial and test functions. Many examples like the one below can be found
throughout the tutorial:

.. code:: python

    u = ufl.TrialFunction(domainSpace)
    v = ufl.TestFunction(rangeSpace)
    op = dune.fem.operator.galerkin(u**2/2*v*dx)
    v_h = rangeSpace.interpolate(...)

    # In the following gf is some grid function,
    # e.g., a function over the domain space:
    op(gf,v_h)
    # To linearize firs construct a linear operator
    linOp = op.linear()
    op.jacobian(gf,linOp) # gf is the function around which to linearize

The call to ``dune.fem.operator.galerkin`` could as mentioned above include
the domain and ranges spaces to use.

- `domainSpace`: the discrete domain space :math:`V_h` for the operator.
- `rangeSpace`:  the discrete range space :math:`W_h` for the operator.
- `linear`:      return a operator class that can store the linearization :math:`DL`.
  This contains an underlying (sparse) matrix structure that
  will depend on the storage argument provided to the spaces.
- `__call__`:    evaluate the operator.
- `jacobian`:    perform the linearization constructing the matrix inside
  the provided linear operator. There are two versions for
  this call: ``op.jacobian(gf,linOp)`` and ``op.jacobian(gf,linOp,b_h)``.
  The first version was described above, the second version
  fills in an addition discrete function ``b_h``. This
  provides the negative of the constant term in the Taylor
  expansion :math:`<L[v_h],w_h> = -<b_h,w_h> + <DL[\bar{v_h}]v_h,w_h>`.
  We provide the negative here but that term is often the
  right hand side needed to solve the problem. For example
  if :math:`L[v_h]` is linear then :math:`DL[0]v_h = b_h` is the linear
  system that defines the solution to :math:`L[v_h]=0`.
- `model`:
- `setCommunicate`:
- `setQuadratureOrders`:
- `gridSizeInterior`:

Operator with DBC:
------------------

If the operator was constructed not only from a form but including
``DirichletBC`` classes which are provided as a tuple/list in the first argument
during construction of the operator:

.. code:: python

    u = ufl.TrialFunction(domainSpace)
    v = ufl.TestFunction(rangeSpace)
    op = dune.fem.operator.galerkin( [u**2/2*v*dx, dune.ufl.DirichletBC(rangeSpace,g)] )

The resulting operator will have the same methods as before but both the
``__call__`` and the ``jacobian`` method will include handing of the
dirichlet boundary conditions:

- `__call__`: this method will as before compute :math:`w_h=L[v_h]` but
  will in addition all components in :math:`w_h` associated with the
  Dirichlet boundary will be of the form :math:`w_i = v_i - g`.

  .. note:: a more mathematical explanation is to consider all functionals
            :math:`\lambda` which are associated with the boundaries, i.e.,
            depend on their argument restricted to the boundary:
            :math:`\lambda(g_1)=\lambda(g_2)` if :math:`g_1=g_2` on the
            dirichlet boundary. After :math:`w_h=L[v_h]` the following holds
            :math:`\lambda(w_h) = \lambda(v_h-g) = \lambda(v_h)-\lambda(g)`
            for each of these functionals.

- `jacobain`: the matrix assembled will contain a unit row for each row
  associated with the Dirichlet constraints.

In addition to the methods from the standard operator these operators have
additional methods:

- `setConstraints`:
- `subConstraints`:
- `dirichletBlocks`:

.. index:: Scheme; Scheme API

Scheme without DBC:
-------------------

In this project we use the term ``scheme`` to describe a special type of
operators with the main additional functionality of solving the system
:math:`L[v_h]=0` (or :math:`L[v_h]=b_h`). A consequence of this is that
the domain and the range space need to be identical - Petrov Galerkin
method which only require that the dimensions match is not currently
available.

All methods from operators are available on schemes. In addition we have

- `space`: this is a convenience method made available since the domain and
  the range space are identical.
- `solve`: this is the central new method available in ``schemes``. 
  A call ``scheme.solve(target=v_h)`` returns a solution to the problem
  :math:`L[v_h]=0`. This is by default based on a Newton-Krylov solver.
- `parameters`:
- `parameterHelp`:
- `preconditioning`:
- `setErrorMeasure`:
- `inverseLinearOperator`: not sure if this is used

Scheme with DBC:
----------------

If a ``scheme`` is constructed including boundary conditions then the
additional method available on the operator are available on the scheme as
well.

Additional
----------

- `virtualization`:

Unclear
-------

- `dimRange`: is available on schemes.
  The following should be removed or on operator (same as ``scheme.rangeSpace.dimRange``
  or even ``scheme.space.dimRange=scheme.domainSpace.dimRange``)
- `gridSizeInterior`: is missing on scheme


