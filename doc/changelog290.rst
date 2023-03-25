Since 2.9 release
=================

Assembly function
-----------------

Improved Pickling support
-------------------------
Discrete functions can now be directly pickled and will
be reconstructed with the grid on load.
The required JIT modules will also be
generated on load if they do not exist on the system.
For this to work one needs to use the ``dune.common.pickle``
module as discussed in the Section on
:ref:`pickling</backuprestore_nb.ipynb#Pickling>`.

Paraview Reader
---------------
A paraview reader is added to the ``dune.fem`` package which can directly work
with the pickled files. UFL transformations
can be applied after loading the discrete functions in paraview.
For paraview to find the reader a environment variable needs to be set.
The correct command can be found running ``python -m dune.fem reader`` or
the correct env variable can be set by running

.. code-block::

  export PV_PLUGIN_PATH=`python -m dune.fem readerpath`

See https://gitlab.dune-project.org/dune-fem/dune-fem/-/merge_requests/581.


Sampler along Boundary
----------------------
In addition to a ``lineSample`` and ``pointSample`` the ``dune.fem.utility``
module now also provides a ``boundarySample`` function

.. code-block::

    boundarySample(gridFunction, boundaryId=None, boundaryDomain=None, order=-1)

See the help text of `dune.fem.utility.boundarySample` for details.

Point sampler returning (entity,local coordinate) for given global coordinate
------------------------------------------------------------------------------

Improved Dirichlet Boundary Indexing in UFL Forms
-------------------------------------------------
To restrict a domain integral to a part of the boundary given by a boundary id 
one can now use ``ufl.ds(id)`` or multiply the integral with ``dune.ufl.BoundaryId(id)``
(see the section on
:ref:`Mixing different types of boundary conditions</boundary_nb.ipynb#Mixing-different-types-of-boundary-conditions>`.

**Note**: we have also unified the
:ref:`default boundary ids</boundary_nb.ipynb#Cartesian-boundary-ids>`
when using a
``dune.grid.cartesianDomain`` to always be the same as when using
``dune.grid.structuredGrid``.
This was not the case with some grid implementations.

Eisenstatwalker in Newton method
--------------------------------

Storage setup in scheme
-----------------------

Methods on space to return local mapper and for evaluating basis functions
--------------------------------------------------------------------------
