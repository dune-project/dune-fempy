# %% [markdown]
# # Checkointing: backup and restore
# We show here how the state of a simulation can be stored and restored
# using pickle. This is still work in progress:
# We implement a class that stores a hierarchical grid and some discrete
# functions and other information.

# %%
import matplotlib
matplotlib.rc( 'image', cmap='jet' )

import pickle, numpy
class CheckPointer:
    # if 'hgrid=None' then restore hierarchical grid using given file name
    # otherwise backup 'hgrid' to given file
    def __init__(self,fileName,hgrid=None):
        if hgrid is None: # restore file
            self.file = open(fileName,'rb')
            self.hgrid = pickle.load(self.file)
        else:             # store file
            self.file = open(fileName,'wb')
            self.hgrid = hgrid
            pickle.dump(self.hgrid,self.file)
        self.items = []

    def hierarchicalGrid(self):
        return self.hgrid

    # add a item to picke when backup is called (or restored)
    def add(self,item):
        self.items += [item]

    def backup(self):
        for i in self.items:
            if hasattr(i,"read"):
                pickle.dump(i.write(),self.file)
            else:
                pickle.dump(i,self.file)
    def restore(self):
        for i in self.items:
            if hasattr(i,"read"):
                i.read( pickle.load(self.file) )
            else:
                i = pickle.load( self.file ) # this doesn't work

# %% [markdown]
# Construct a grid, do some local adaptivity based on gradient of discrete
# function, and store to file.

# %%
import math
from dune.alugrid import aluConformGrid as view
from dune.grid import gridFunction, cartesianDomain, Marker
from dune.fem.space import lagrange
from dune.fem import adapt

def markh(element):
    center = element.geometry.referenceElement.center
    eta    = estimate.localFunction(element).evaluate(center)[0]
    if eta > hTol and element.level < maxLevel:
      return Marker.refine
    elif eta < 0.01*hTol:
      return Marker.coarsen
    else:
      return Marker.keep

def test1(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[10,10]) )

    @gridFunction(grid)
    def gf(x): return math.sqrt(1-(x[0]**2+x[1]**2)) if x[0]**2+x[1]**2<1 else 0

    space = lagrange(grid)
    df = space.interpolate(gf,name="test")

    grid.hierarchicalGrid.globalRefine(4)
    for i in range(5):
        grid.hierarchicalGrid.mark(lambda e:
             Marker.refine if df.localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 1
             else Marker.coarsen)
        adapt(grid.hierarchicalGrid)
        df.interpolate( gf )
    print("size of adapted grid:", grid.size(0))
    df.plot()

    checkPointer = CheckPointer(fileName, grid.hierarchicalGrid)
    checkPointer.add(df)
    checkPointer.backup()

test1("dump.pickle")

# %% [markdown]
# Now we restore the backup

# %%
checkPointer = CheckPointer("dump.pickle")
grid = checkPointer.hierarchicalGrid().leafView
print("size of restored grid:", grid.size(0))
space = lagrange(grid)
df = space.interpolate([0],name="test")
checkPointer.add(df)

checkPointer.restore()
df.plot()
