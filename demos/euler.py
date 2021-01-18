"""
Created on Thu Jun 11 14:52:41 2020

modified from dune-fempy demo:
https://www.dune-project.org/sphinx/content/sphinx/dune-fem/euler_nb.html (May 2020)
@author: Peter Meisrimel, Lund University, with lots of help from Robert Klöfkorn
"""

## Euler System of Gas Dynamics

import numpy as np
import ufl
from dune.grid import reader
from dune.fem.space import dgonb
from dune.femdg import femDGOperator
from dune.femdg.rk import femdgStepper
from dune.alugrid import aluSimplexGrid
from dune.ufl import expression2GF, Constant
from dune.fem.utility import Sampler
from scipy.interpolate import interp1d

class Model:
    def __init__(self):
        ## material parameters
        ## source https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
        self.cp, self.cv = 1005, 718
        self.gamma = 1.4
        self.R = 287.058 ## R = cp - cv ## specific gas constant

        ## initial conditions
        self.rho0 = 1.225 # kg/m**3
        self.T0 = 293.15 # 20 C
        self.T1 = 900
        self.vx0, self.vy0 = 1, 0 # m/s

        ## initial energies
        self.E0 = self.rho0 * (self.R*self.T0*(self.gamma - 1) + 0.5 * (self.vx0**2 + self.vy0**2))
        self.E1 = self.rho0 * (self.R*self.T1*(self.gamma - 1) + 0.5 * (self.vx0**2 + self.vy0**2))

        ## by ID, see triangle.dgf file, 5 = interface
        # given in conservative form
        self.boundary = {1: lambda t, x, U: U, ## right
                         2: lambda t, x, U: [self.rho0, self.vx0, self.vy0, self.E0], ## upper
                         3: lambda t, x, U: [self.rho0, self.vx0, self.vy0, self.E0], ## left
                         4: lambda t, x, U: [U[0], U[1], 0,
                                             U[0]*(self.R*(self.gamma-1)*self.T0 + 0.5*(U[1]/U[0])**2)], ## lower left
                         5: lambda t, x, U: [U[0], U[1], 0,
                                             U[0]*(self.R*(self.gamma-1)*self.T1 + 0.5*(U[1]/U[0])**2)], ## lower centre
                         6: lambda t, x, U: [U[0], U[1], 0,
                                             U[0]*(self.R*(self.gamma-1)*self.T0 + 0.5*(U[1]/U[0])**2)]  ## lower right
                         }

        self.eps = 1e-13
        self.x = ufl.SpatialCoordinate(ufl.triangle)
        # ufl.conditional(cond, True, False)
        self.u0 = ufl.conditional(ufl.le(self.x[1], -1 + self.eps), ## condition for lower boundary
                          ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E1]),
                          ufl.as_vector([self.rho0, self.vx0, self.vy0, self.E0]))

        self.n_q = Constant((0., -1.))

        ## p. 259 http://inis.jinr.ru/sl/Simulation/Tannehill,_CFM_and_Heat_Transfer,2_ed/chap05.pdf
        self.C3 = 2.495 * 10**(-3)
        self.C4 = 194

    # helper function
    def toPrim(self, U):
        v = ufl.as_vector([U[i]/U[0] for i in range(1,3)])
        kin = ufl.dot(v,v) * U[0]/2
        pressure = (self.gamma - 1)*(U[3] - kin)
        return U[0], v, pressure

    # interface methods for model
    def F_c(self, t, x, U):
        rho, v, p = self.toPrim(U)
        return ufl.as_matrix([
                  [rho*v[0], rho*v[1]],
                  [rho*v[0]*v[0] + p, rho*v[0]*v[1]],
                  [rho*v[0]*v[1], rho*v[1]*v[1] + p],
                  [(U[3] + p)*v[0], (U[3] + p)*v[1]]])

    # interface method needed for LLF and time step control
    def maxLambda(self, t, x, U, n):
        rho, v, p = self.toPrim(U)
        return abs(ufl.dot(v,n)) + ufl.sqrt(self.gamma*p/rho)

    ## methods for limiting
    def velocity(self, t, x, U):
        _, v, _ = self.toPrim(U)
        return v
    def physical(self, t, x, U):
        rho, _, p = self.toPrim(U)
        return ufl.conditional(rho > 1e-8, ufl.conditional(p > 1e-8, 1, 0), 0)
    def jump(self, t, x, U, V):
        _, _, pL = self.toPrim(U)
        _, _, pR = self.toPrim(V)
        return (pL - pR)/(0.5*(pL + pR))

    def temp(self, U):
        rho, _, p = self.toPrim(U)
        return p/(rho*self.R)

    def thermal_cond(self, U):
        T = self.temp(U)
        return self.C3 * T**(1.5)/(T + self.C4)

    def set_temp(self, T, U):
        rho, v, p = self.toPrim(U)
        return ufl.as_vector((U[0], U[1], 0, rho*(self.R*(self.gamma - 1)*T + 0.5*v[0]*v[0])))

## TODO
"""
Change domain/grid
preconditioning, is it already included?
upgrade to navier stokes?
add functionality such that an initial step is done during initialization to build necessary dune functions
possibly compute steady state solution by running until difference to previous timestep is smaller than a given tolerances,
rather than running a fixed number of steps
"""

class Problem_heat_euler_fluid:
    u0_filename = "u0.out"
    ## get some more inputs, e.g. various orders
    def __init__(self, num_points_interface = 32, order_space = 2, order_time = 2, load_u0_from_file = True):
        self.Model = Model()

        self.domain = (reader.dgf, "triangle.dgf") ## read domain specifications from file
        self.gridView = aluSimplexGrid(self.domain, serial = True, dimgrid = 2)
        self.space = dgonb(self.gridView, dimRange = 4, storage = 'istl')
        self.uh = self.space.interpolate(self.Model.u0, name = "solution")
        if load_u0_from_file:
            self.load_u0_from_file()

        self.u_checkpoint = self.space.interpolate(self.Model.u0, name = "checkpoint")

        self.uh_boundary = self.space.interpolate((self.Model.rho0, self.Model.vx0, self.Model.vy0, self.Model.E1), name = "uh_bd")
        self.Model.boundary[5] = lambda t, x, U: self.uh_boundary

        #############
        from dune.fem import parameter
        parameter.append({"fem.verboserank": 0})

        para = {"fem.ode.odesolver": "EX",
                "fem.ode.order": 1,
                "fem.ode.verbose": "full",
                "fem.solver.newton.linear.verbose": True,
                "fem.solver.newton.verbose": True,
#                "fem.solver.newton.linear.preconditioning.method": "ilu"
                }
        self.operator = femDGOperator(self.Model, self.space) # prev: MinMod
        self.stepper = femdgStepper(operator = self.operator, parameters = para)
        self.operator.applyLimiter(self.uh)

        ##############

        ## flux computation
        self.n_q = Constant((0., -1.))
        self.heat_flux = expression2GF(self.uh.space.grid,
                                       self.Model.thermal_cond(self.uh)*ufl.dot(ufl.grad(self.Model.temp(self.uh)), self.n_q),
                                       self.uh.space.order)

        self.NN = num_points_interface + 2 # number of sampling points
        self.sample_points = np.linspace(-0.5, 0.5, self.NN) ## bottom boundary line

        self.sampler = Sampler(self.heat_flux)
        self.get_flux_sample = lambda : self.sampler.lineSample([-0.5, -1], [0.5, -1], self.NN)[1]

        self.intf_rho_eval = expression2GF(self.uh.space.grid, ufl.dot(self.uh, Constant((1, 0, 0, 0))), self.uh.space.order)
        self.rho_sampler = Sampler(self.intf_rho_eval)
        self.get_rho_sample = lambda : self.rho_sampler.lineSample([-0.5, -1], [0.5, -1], self.NN)[1]

        self.intf_rho_vx_eval = expression2GF(self.uh.space.grid, ufl.dot(self.uh, Constant((0, 1, 0, 0))), self.uh.space.order)
        self.rho_vx_sampler = Sampler(self.intf_rho_vx_eval)
        self.get_rho_vx_sample = lambda : self.rho_vx_sampler.lineSample([-0.5, -1], [0.5, -1], self.NN)[1]

        self.temp_fac = Constant(self.Model.R*(self.Model.gamma - 1))
        self.t = 0. ## TODO: possibly needs better treatment for creating checkpoint and resetting

    def load_u0_from_file(self):
        import dune.generator.algorithm as algorithm
        try:
            readData = algorithm.load("readData", "dataio.hh", self.uh, self.u0_filename)
        except:
            readData = algorithm.load("readData", "dune_python/dataio.hh", self.uh, self.u0_filename)
        readData(self.uh, self.u0_filename)

    def write_u0_to_file(self):
        import dune.generator.algorithm as algorithm
        try:
            writeData = algorithm.load("writeData", "dataio.hh", self.uh, self.u0_filename)
        except:
            writeData = algorithm.load("writeData", "dune_python/dataio.hh", self.uh, self.u0_filename)
        writeData(self.uh, self.u0_filename)

    def create_checkpoint(self):
        ## TODO: won't work with macrostepping as of now. Requires storing and loading, possibly via file?
        self.u_checkpoint.interpolate(self.uh)

    def reset(self):
        self.uh.interpolate(self.u_checkpoint)
        self.t = 0.

    def get_u0(self):
        return self.get_flux()

    def get_flux(self):
        return self.get_flux_sample()

    def set_intf_temp(self, temps):
        intp_temp = interp1d(self.sample_points, temps, fill_value = 'extrapolate')
        intp_rho = interp1d(self.sample_points, self.get_rho_sample(), fill_value = 'extrapolate')
        intp_rho_vx = interp1d(self.sample_points, self.get_rho_vx_sample(), fill_value = 'extrapolate')

        def f(x):
            rho = intp_rho(x)
            vx = intp_rho_vx(x)
            return (rho, vx, 0, rho*(self.temp_fac*intp_temp(x) + 0.5*vx**2/rho))

        self.uh_boundary.interpolate(lambda x: f(x[0]))

    def do_step(self, t, dt, vals):
        self.operator.setTime(t)

        self.set_intf_temp(vals)
        self.t += self.stepper(self.uh, dt = dt)

        return self.get_flux()

    def solve(self, tf, n):
        dt = tf/n
        self.t = 0
        for i in range(n):
            print('doing step', self.t, dt, flush = True)
            vals = np.array([900.]*self.NN)
            self.do_step(self.t, dt, vals)
            ## too large timesteps lead to arbortian of solving process

            from dune.grid import OutputType
            self.gridView.writeVTK('solution', pointdata={'T': self.Model.temp(self.uh)}, number = i, outputType=OutputType.ascii)

if __name__ == '__main__':
    tf = 0.05
    n = 100
    solver = Problem_heat_euler_fluid(load_u0_from_file = False)
    solver.create_checkpoint()
    solver.solve(tf, n)
    solver.write_u0_to_file()
