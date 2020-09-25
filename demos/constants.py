# <markdowncell>

# TODO document the problem here

# <codecell>

import os
# ensure some compilation output for this example
os.environ['DUNE_LOG_LEVEL'] = 'info'
print("Using DUNE_LOG_LEVEL=",os.getenv('DUNE_LOG_LEVEL'))

from dune.grid import structuredGrid as leafGridView
from dune.ufl import Constant
from ufl import SpatialCoordinate, triangle
from dune.fem.function import integrate, uflFunction

import time

gridView = leafGridView([0, 0], [1, 1], [4, 4])

x = SpatialCoordinate(triangle)

def computeArea( function, t ):
    xarea = integrate(gridView, function, order=1)
    err = abs(xarea - float(t)*0.5)
    if err > 1e-5:
        print("Wrong result for integral computation, error = ",err)
    print("t=",float(t), " integral = ", xarea )

dt = 1.

#############################################################
#   Example 1: The naive way
#############################################################
# use system time as initializer to force re-compilation
# because of the changing time
t = time.time()
exact = t*x[0]

endTime = t+5.
# take start time for later measurement
startTime = time.time()

# Wrong because t in exact is always initial t since
# that is "hard-coded" into the ufl description
while t<endTime:
    computeArea( exact, t )
    t += dt

print("Finished, CPU time in sec: ",time.time() - startTime)

#############################################################
#   Example 2: Use dune.ufl.Constant, but no uflFunction
#############################################################

# Make t a dune.ufl.Constant to be able to update the value
# use t.value = ... to change the actual value
t = Constant(time.time(), name="t")
exact = t*x[0]
endTime = t.value+5.

# take start time for later measurement
startTime = time.time()

# Now the result will be right, but the change of the value of t
# will cause re-compilations in every step
while t.value<endTime:
    xarea = integrate(gridView, exact, order=1)
    print("t=",t.value, " integral = ", xarea )
    t.value += dt

# before output, make sure that all quantities are computed even after grid chan
print("Finished, CPU time in sec: ",time.time() - startTime)

#############################################################
#   Example 3: Use dune.ufl.Constant and uflFunction
#############################################################

# Again use a dune.ufl.Constant to be able to update the value
# use t.value = ... to change the actual value
t = Constant(time.time(), name="t")
exact = t*x[0]
# convert ufl expression into a ufl gridFunction, by doing this
# Constants are introduced as dynamically changable variables
gf = uflFunction(gridView, name="exact", order=1, ufl=exact)
endTime = t.value+5.
# take start time for later measurement
startTime = time.time()
while t.value<endTime:
    computeArea( gf, t )
    t.value += dt

# before output, make sure that all quantities are computed even after grid chan
print("Finished, CPU time in sec: ",time.time() - startTime)
