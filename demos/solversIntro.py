import numpy, sys, io
import matplotlib.pyplot as plt

from dune.grid import structuredGrid as leafGridView
from dune.fem.space import lagrange
from dune.fem.function import gridFunction
from dune.fem import integrate
from dune.ufl import Constant
from ufl import TestFunction, TrialFunction, SpatialCoordinate, FacetNormal, \
                dx, ds, div, grad, dot, inner, sqrt, exp, sin,\
                conditional

gridView = leafGridView([0, 0], [1, 1], [12, 12])
gridView.hierarchicalGrid.globalRefine(1)
space = lagrange(gridView, order=2)

x = SpatialCoordinate(space)
initial = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
exact   = lambda t: exp(-2*t)*(initial - 1) + 1

u_h   = space.interpolate(initial, name='u_h')
u_h_n = u_h.copy(name="previous")

u = TrialFunction(space)
v = TestFunction(space)
dt = Constant(0.025, name="dt")   # time step
t  = Constant(0, name="time")     # current time

abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2/(1 + sqrt(1 + 4*abs_du(u)))
a = ( dot((u - u_h_n)/dt, v) \
    + 0.5*dot(K(u)*grad(u), grad(v)) \
    + 0.5*dot(K(u_h_n)*grad(u_h_n), grad(v)) ) * dx

f = lambda s: -2*exp(-2*s)*(initial - 1) - div( K(exact(s))*grad(exact(s)) )
g = lambda s: K(exact(s))*grad(exact(s))
n = FacetNormal(space)
bf = 0.5*(f(t)+f(t+dt))*v*dx
bg = 0.5*dot(g(t)+g(t+dt),n)*v*ds

import time as timer
def evolve(scheme, u_h, u_h_n, endTime):
    start = timer.time()
    time = 0
    nIter,linIter = 0,0
    steps = 0
    while time < (endTime - 1e-6):
        scheme.model.time = time
        t.value = time
        # assert t.value == time
        u_h_n.assign(u_h)
        info = scheme.solve(target=u_h)
        nIter += info["iterations"]
        linIter += info["linear_iterations"]
        time += scheme.model.dt
        steps += 1
    print(f"steps={[steps,scheme.model.dt]}, runtime={timer.time()- start}, iterations={[nIter,linIter]}",flush=True)

endTime    = 0.25
exact_end  = exact(endTime)
l2error = gridFunction(name="l2error", expr=dot(u_h - exact_end, u_h - exact_end))
h1error = gridFunction(name="h1error", expr=dot(grad(u_h - exact_end), grad(u_h - exact_end)))
