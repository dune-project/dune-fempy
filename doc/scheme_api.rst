API for schemes and operators
=============================

Operator without DBC:
---------------------

- domainSpace:
- gridSizeInterior:
- jacobian:
- linear:
- model:
- rangeSpace:
- setCommunicate:
- setQuadratureOrders:

Operator with DBC:
------------------

- dirichletBlocks:
- setConstraints:
- subConstraints:

Scheme without DBC:
-------------------

- From operator:
   - domainSpace:
   - jacobian:
   - linear:
   - model:
   - rangeSpace:
   - setCommunicate:
   - setQuadratureOrders:
- Should be removed or on operator (same as scheme.rangeSpace.dimRange or even scheme.space.dimRange=scheme.domainSpace.dimRange)
   - dimRange:
- Missing on scheme:
   - gridSizeInterior:
- New on scheme:  
   - inverseLinearOperator:
   - parameterHelp:
   - parameters:
   - preconditioning:
   - setErrorMeasure:
   - solve:
   - space: (since domainSpace=rangeSpace)

Scheme with DBC:
----------------

- From operator with DBC:
   - dirichletBlocks:
   - setConstraints:
   - subConstraints:
