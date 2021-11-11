.. index:: Grid construction; Filtered grid views

FilteredGridView examples
=========================

In this example we use ``filteredGridView`` to filter out parts of a
given grid view.

.. code:: ipython3

    
    from ufl import SpatialCoordinate, dot
    from dune.grid import cartesianDomain
    from dune.alugrid import aluConformGrid as leafGridView
    from dune.fem.view import filteredGridView
    from dune.fem.space import lagrange

Create a host grid view of the underlying grid as usual.

.. code:: ipython3

    
    gridView = leafGridView( cartesianDomain([0,0],[1,1],[16,16]) )


.. parsed-literal::

    
    Created parallel ALUGrid<2,2,simplex,conforming> from input stream. 
    


.. parsed-literal::

    GridParameterBlock: Parameter 'refinementedge' not specified, defaulting to 'ARBITRARY'.
    WARNING (ignored): Could not open file 'alugrid.cfg', using default values 0 < [balance] < 1.2, partitioning method 'ALUGRID_SpaceFillingCurve(9)'.
    
    You are using DUNE-ALUGrid, please don't forget to cite the paper:
    Alkaemper, Dedner, Kloefkorn, Nolte. The DUNE-ALUGrid Module, 2016.


Now create a ``filteredGridView`` with a simple callable specifying
which element belongs to the new domain.

.. code:: ipython3

    
    filteredView = filteredGridView(gridView, lambda e: e.geometry.center.two_norm > 0.5, domainId=1)
    space = lagrange(filteredView, order=2)
    x = SpatialCoordinate(space)
    solution = space.interpolate(dot(x,x),name="solution")
    solution.plot()
    print("number of dofs:", solution.size,\
          "integral over filtered domain",solution.integrate())



.. image:: filteredgridview_nb_files/filteredgridview_nb_5_0.png


.. parsed-literal::

    number of dofs: 1089 integral over filtered domain 0.6413014729817705


Revert the filter and create a ``filteredGridView`` with also provides
an overloaded index set with a consecutive index for the entities
belonging to the new domain.

.. code:: ipython3

    filteredView = filteredGridView(gridView, lambda e: e.geometry.center.two_norm < 0.5, domainId=1,
                                    useFilteredIndexSet=True)
    space = lagrange(filteredView, order=2)
    x = SpatialCoordinate(space)
    solution = space.interpolate(dot(x,x),name="solution")
    solution.plot()
    print("number of dofs:", solution.size,\
          "integral over filtered domain",solution.integrate())



.. image:: filteredgridview_nb_files/filteredgridview_nb_7_0.png


.. parsed-literal::

    number of dofs: 239 integral over filtered domain 0.025365193684895825


Compare with the original grid view.

.. code:: ipython3

    space = lagrange(gridView, order=2)
    solution = space.interpolate(dot(x,x),name="solution")
    solution.plot()
    print("number of dofs:", solution.size,\
          "integral over filtered domain",solution.integrate())



.. image:: filteredgridview_nb_files/filteredgridview_nb_9_0.png


.. parsed-literal::

    number of dofs: 1089 integral over filtered domain 0.6666666666666667

