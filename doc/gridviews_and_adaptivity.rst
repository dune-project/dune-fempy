.. sectionauthor:: Andreas Dedner <a.s.dedner@warwick.ac.uk>, Robert Kl\ |oe|\ fkorn

.. index:: Adaptation; Adaptation overview
.. index:: GridView; GridView examples

###############################################
Grid views and adaptivity
###############################################

=======================================================
Overview and some basic grid views (level and filtered)
=======================================================

When constructing a grid, the object returned to Python is always the so
called ``LeafGridView``. Without any refinement this is simply a view on all
the elements of the grid. As soon as the grid is refined the leaf grid view changes
so that it always contains the ``leaf`` elements the grid, i.e., the elements
on the finest level. Since it is a read only view refinement is carried out
using the underlying ``hierarchical grid``, i.e.,

.. code:: python

   grid.hierarchicalGrid.globalRefine(1)

For a given hierarchical grid one can use different views, i.e., a view on
all the elements of a given level:

.. code:: python

   levelView = grid.hierarchicalGrid.levelView(1)

.. toctree::
   :maxdepth: 2
   :name: levelGV

   levelgridview_nb


DUNE-FEM provides a number of additional views which will be discussed
further in this chapter:

*  ``dune.fem.view.adaptiveLeafGridView``: this view should be used when
   the grid is supposed to be locally adapted. The view is still on the leaf
   elements of the grid but way data is attached to the entities of the grid
   is optimized for frequent grid changes as caused by local adaptivity. Its
   usage is shown in the following examples.

*  ``dune.fem.view.geometryGridView``: this is an example of a ``meta`` grid
   view which is constructed from an existing grid view and replaces some
   aspect - in this case the geometry of each element using a given grid
   function. This concept makes it easy to perform simulations for example on
   complex domains or on moving grids as shown in :ref:`Evolving Domains<geomGV>`.

*  ``dune.fem.view.filteredGridView``: allows the user to select a subset of
   a given grid view by providing a filter on the elements. In its simplest
   version only the iterator is replaced with an iterator over the elements
   in the filter but in addition it is also possible to obtain a new index
   set with indices restricted to the elements in the filter.

.. toctree::
   :maxdepth: 2
   :name: filteredGV

   filteredgridview_nb

===========================================================
Dynamic Local Grid Refinement and Coarsening (h-adaptation)
===========================================================

For refining and coarsening a grid locally the ``dune.fem`` module provides a
function ``gridAdapt``. The storage of **all** discrete functions will be
automatically resized to accommodate the changes in the hierarchical grid those
are associated with.

.. attention::
    The resulting DOF vector will **not be initialized**. To prolong and restrict data
    from the old to the new grid, the corresponding discrete functions have to be
    passed to the ``dune.fem.gridAdapt`` method:


.. code:: python

    fem.gridAdapt(marker, [u1,u2,...,uN]) # marker is explained below

.. note::

    All the prolongation/restriction described here **requires** that the
    discrete spaces were constructed over an ``adaptiveLeafGridView``.

So the code has to be modified for example as follows

.. code:: python

    from dune.alugrid import aluConformGrid as leafGridView
    from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
    gridView = adaptiveGridView( leafGridView(domain) )

.. attention::

   if the underlying storage of a discrete function is stored on the Python
   side as a numpy array, i.e., ``vec = uh.as_numpy`` was called, then access
   to ``vec`` will be undefined after a grid modification since the
   underlying buffer change will not have been registered.


.. note::
    The marking of elements for refinement or coarsening is provided the by ``marker`` which
    is either an instance of ``GridMarker`` or a callback (see below).

The module ``dune.fem`` provides an object for marking elements for
refinement/coarsening:

.. code:: python

    from dune.fem import GridMarker
    marker = GridMarker( indicator, # grid function to be evaluated at cell center
                         refineTolerance, coarsenTolerance=0, # tolerances (float or callable returning float)
                         minLevel=0, maxLevel=None,    # min and max level allowed
                         minVolume=-1., maxVolume=-1., # min and max element volume allowed
                         statistics = False,           # element counts for marked refined or coarsen
                         strategy = 'default'          # marking strategy: default or doerfler
                         layered = 0.05,               # parameter for Doerfler layered strategy.
                         markNeighbors = False )       # for elements marked for refinement
                                                       # also mark neighboring elements

where ``indicator`` is a piecewise constant grid function.
An element :math:`E` is marked for refinement if the value of ``indicator``
evaluated at the center of :math:`E` is greater than ``refineTolerance`` and coarsened if the
value is less than ``coarsenTolerance``. The element :math:`E` is not
refined if its level is already at ``maxLevel`` and not coarsened if its
level it at ``minLevel`` and accordingly for the ``minVolume`` and ``maxVolume``.
This method can for example be used to refine the grid according to an equal distribution strategy by
setting ``refineTolerance=theta/grid.size(0)`` with some tolerance ``theta``.
A **layered Doerfler strategy** is also available. Select ``strategy='doerfler'``
in the above ``GridMarker`` object constructor.

.. note::
   For debugging and testing the marking of elements can also be done by
   providing a callback to the ``gridAdapt`` funtion with the following
   signature:

.. code:: python

    from dune.grid import Marker
    def marker( element ) -> Marker
        if refine element: # if condition for refinement true
            return Marker.Refine
        elif coarsen element: # if condition for coarsening true
            return Marker.coarsen
        else: # otherwise keep element as is (might still be changed)
            return Marker.keep

The module ``dune.fem`` also provides a ``globalRefine(level,*dfs)`` method,
where a negative level globally coarsens the grid. If discrete functions
are passed in they will be prolonged (restricted) and resized correctly,
the DOF vectors of all other discrete functions will only be resized, i.e.
stored data will be lost.

The following two examples showcase adaptivity: the first one using a
residual a-posteriori estimator for an elliptic problem, the second one
shows adaptivity for a time dependent phase field model for crystal growth.
At the end of this section a dual weighted residual approach is used to
optimize the grid with respect to the error at a given point. While the
first two examples can be implemented completely using the available Python
bindings the final example requires using a small C++ snippet which is easy
to integrate into the Python code.

.. toctree::
   :maxdepth: 2

   laplace-adaptive_nb
   crystal_nb

=================================================
Overview on adaptivity for spaces (p-adaptation)
=================================================

For adjusting the polynomial degree of a discrete function space locally (on each element separately)
the ``dune.fem`` module provides a function ``spaceAdapt``. The storage of **all** discrete functions created over this space
will be automatically resized to accommodate the changes in the space that is adapted.

.. attention::
    Unlike for grid adaptation the **space adaptation** has to be carried out **for each space separately**.
    To correctly transfer data the corresponding discrete functions have to be
    passed to the ``dune.fem.spaceAdapt`` method:


.. code:: python

    fem.spaceAdapt(marker, [u1,u2,...,uN]) # marker is explained below

.. note::

    All the prolongation/restriction described here **requires** that the
    discrete spaces were constructed over an ``adaptiveLeafGridView``.

    Furthermore, the discrete space used has to **support p-adaptation**, which is
    only the case for a handfull of space. This can be checked with
    ``space.canAdapt``, which is true if the space supports p-adaptation.

So the code has to be modified for example as follows

.. code:: python

    from dune.alugrid import aluConformGrid as leafGridView
    from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
    gridView = adaptiveGridView( leafGridView(domain) )

.. attention::

   if the underlying storage of a discrete function is stored on the Python
   side as a numpy array, i.e., ``vec = uh.as_numpy`` was called, then access
   to ``vec`` will be undefined after a grid modification since the
   underlying buffer change will not have been registered.


.. note::
    The marking of the local polynomial degree of a discrete space is provided the by ``pmarker`` which
    is either an instance of ``SpaceMarker`` or a callback (see below).

The module ``dune.fem`` provides an object for marking spaces for
p-refinement/p-coarsening:

.. code:: python

    from dune.fem import SpaceMarker
    pmarker = SpaceMarker( indicator, # grid function to be evaluated at cell center
                           refineTolerance, coarsenTolerance=0, # tolerances (float or callable returning float)
                           minOrder=0,  # minimal order to be assumed (in [0, space.order)
                           maxOrder=-1, # maximal order to be assumed (in [minOrder, space.order])
                           markNeighbors = False, # for elements marked for increase also mark neighboring elements
                           statistics = False)    # element counts for marked refined or coarsen

where ``indicator`` is a piecewise constant grid function.
The polynomial degree :math:`p` of a space for an element :math:`E` is increased by :math:`1`
if the value of ``indicator`` evaluated at the center of :math:`E` is greater than ``refineTolerance`` and decreased by :math:`1`
if the value is less than ``coarsenTolerance``. The polynomial degree remains of
a space for :math:`E` is not changed if the result would be outside of :math:`[minOrder, maxOrder]`.

.. note::
   For debugging and testing the marking of elements can also be done by
   providing a callback to the ``spaceAdapt`` funtion with the following
   signature:

.. code:: python

    def pmarker( element ) -> int
        # current polynomial degree
        p = space.localOrder(element)
        if increase degree:
            return p+1
        elif decrease degree: # if condition for coarsening true
            return p-1
        else: # otherwise keep current degree
            return p

An example where p-adaptation is used can be found here:

.. toctree::
   :maxdepth: 2

   twophaseflow_descr


=================================================
Overview on moving grids (r-adaptation)
=================================================

An example of `r-adaptation` for a two-phase flow problem in porous media is found
on the `Dune-MMesh`_ documentation page.

.. _Dune-MMesh: https://dune-mmesh.readthedocs.io/en/latest/examples/moving.html

As mentioned above DUNE-FEM provides a grid view that makes it easy to
exchange the geometry of each entity in the grid. To setup such a grid view
one first needs to construct a standard grid view, i.e., a ``leafGridView``
and define a grid function over this view using for example a discrete
function, a UFL function, or one of the concepts described in the section
:ref:`Grid Function</concepts_nb.ipynb#Grid-Functions>`.
Note that the topology of the
grid does not change, i.e., how entities are connected with each other.
The following shows an example of how to change a grid of the unit square
into a grid of a diamond shape:

.. include:: geoview_nb.rst

By using a discrete function to construct a geometry grid view, it becomes
possible to simulate problems on evolving domains where the evolution is
itself the solution of the partial differential equation. We demonstrate
this approach based on the example of surface mean curvature flow first in
its simplest setting and then with the evolution of the surface depending
on values of a computed surface quantity satisfying a heat equation on the
surface:

Further examples using moving meshes:


.. toctree::
   :maxdepth: 2
   :name: geomGV

   mcf_nb
   elasticity_nb


