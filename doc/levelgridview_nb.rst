.. code:: ipython3

    from dune.grid import cartesianDomain
    from dune.alugrid import aluConformGrid, aluSimplexGrid
    domain = cartesianDomain([0,0],[1,1],[1,1])

.. code:: ipython3

    # first construct a grid using quartering as refinement strategy
    aluView = aluSimplexGrid(domain)
    hGrid = aluView.hierarchicalGrid
    hGrid.globalRefine(5)
    for level in range(hGrid.maxLevel):
        print("level:",level, "number of elements:",hGrid.levelView(level).size(0))


.. parsed-literal::

    
    Created parallel ALUGrid<2,2,simplex,nonconforming> from input stream. 
    


.. parsed-literal::

    WARNING (ignored): Could not open file 'alugrid.cfg', using default values 0 < [balance] < 1.2, partitioning method 'ALUGRID_SpaceFillingCurve(9)'.
    
    You are using DUNE-ALUGrid, please don't forget to cite the paper:
    Alkaemper, Dedner, Kloefkorn, Nolte. The DUNE-ALUGrid Module, 2016.


.. parsed-literal::

    level: 0 number of elements: 2
    level: 1 number of elements: 8
    level: 2 number of elements: 32
    level: 3 number of elements: 128
    level: 4 number of elements: 512


.. code:: ipython3

    # now construct a grid using bisection as refinement strategy
    aluView = aluConformGrid(domain)
    hGrid = aluView.hierarchicalGrid
    hGrid.globalRefine(10)
    for level in range(hGrid.maxLevel):
        print("level:",level, "number of elements:",hGrid.levelView(level).size(0))


.. parsed-literal::

    GridParameterBlock: Parameter 'refinementedge' not specified, defaulting to 'ARBITRARY'.
    WARNING (ignored): Could not open file 'alugrid.cfg', using default values 0 < [balance] < 1.2, partitioning method 'ALUGRID_SpaceFillingCurve(9)'.


.. parsed-literal::

    
    Created parallel ALUGrid<2,2,simplex,conforming> from input stream. 
    
    level: 0 number of elements: 2
    level: 1 number of elements: 4
    level: 2 number of elements: 8
    level: 3 number of elements: 16
    level: 4 number of elements: 32
    level: 5 number of elements: 64
    level: 6 number of elements: 128
    level: 7 number of elements: 256
    level: 8 number of elements: 512
    level: 9 number of elements: 1024

