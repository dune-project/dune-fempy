from ufl import *

from dune.grid import cartesianDomain, gridFunction
from dune.alugrid import aluConformGrid as leafGridView
from dune.fem.view import adaptiveLeafGridView as adaptiveGridView
from dune.fem.view import filteredGridView
from dune.fem import parameter
from dune.ufl import Space

import dune.create as create

parameter.append({"fem.verboserank": 0})

domain = cartesianDomain([0,0],[1,1],[16,16])
gridView = adaptiveGridView( leafGridView(domain) )

def contains( e ):
    return e.geometry.center[ 0 ]  > 0.5

filteredView = filteredGridView(gridView, contains, domainId=0 )
space = create.space("lagrange", filteredView, order=2)

solution = space.interpolate([1],name="solution")
filteredView.writeVTK("filtered", pointdata=[solution])
