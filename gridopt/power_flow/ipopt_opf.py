#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import numpy as np
from .method_error import *
from .method import PFmethod

class IpoptOPF(PFmethod):
    """
    IPOPT-based optimal power flow method.
    """
 
    name = 'IpoptOPF'

    parameters = {'weight_cost': 1e0,     # for generation cost
                  'weight_mag_reg': 0.,   # for soft limits
                  'weight_ang_reg': 0.,   # for voltage angle regularization
                  'weight_gen_reg': 0.}   # for generators regularization

    def __init__(self):

        from optalg.opt_solver import OptSolverIpopt
        
        PFmethod.__init__(self)
        parameters = OptSolverIpopt.parameters.copy()
        parameters.update(IpoptOPF.parameters)
        self.parameters = parameters

    def create_problem(self,net):

        import pfnet

        # Parameters
        params = self.parameters
        wc = params['weight_cost']
        wm = params['weight_mag_reg']
        wa = params['weight_ang_reg']
        wg = params['weight_gen_reg']
        
        # Clear flags
        net.clear_flags()
        
        # Voltage magnitudes
        net.set_flags('bus',
                      ['variable','bounded'], 
                      'any',
                      'voltage magnitude')
        
        # Voltage angles
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')

        # Generator powers
        net.set_flags('generator',
                      ['variable','bounded'],
                      'not on outage',
                      ['active power','reactive power'])

        try:
            assert(net.num_vars == (2*net.num_buses-net.get_num_slack_buses() +
                                    2*net.get_num_gens_not_on_outage())*net.num_periods)
            assert(net.num_bounded == (2*net.get_num_gens_not_on_outage() +
                                       net.num_buses)*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
                                    
        # Set up problem
        problem = pfnet.Problem(net)
        problem.add_constraint(pfnet.Constraint('AC power balance',net))
        problem.add_constraint(pfnet.Constraint('variable bounds',net)) 
        problem.add_function(pfnet.Function('generation cost',wc/max([net.num_generators,1.]),net))
        if wm:
            problem.add_function(pfnet.Function('soft voltage magnitude limits',wm/max([net.num_buses,1.]),net))
        if wa:
            problem.add_function(pfnet.Function('voltage angle regularization',wa/max([net.num_buses,1.]),net))
        if wg:
            problem.add_function(pfnet.Function('generator powers regularization',wg/max([net.num_generators,1.]),net))
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net):

        from optalg.opt_solver import OptSolverError, OptSolverIpopt
        
        # Parameters
        params = self.parameters

        # Problem
        problem = self.create_problem(net)
                
        # Set up solver
        solver = OptSolverIpopt()
        solver.set_parameters(params)
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:

            # Update network properties
            net.update_properties(solver.get_primal_variables())
            
            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(solver.get_primal_variables())
            self.set_dual_variables(solver.get_dual_variables())
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)

            # Restors net properties
            net.update_properties()

    def update_network(self,net):
        
        # Get data
        problem = self.results['problem']
        x = self.results['primal variables']
        lam,nu,mu,pi = self.results['dual variables']
       
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == x.shape)
        assert(net.num_vars == x.size)
        assert(problem.A.shape[0] == lam.size)
        assert(problem.f.shape[0] == nu.size)
        assert(problem.G.shape[0] == mu.size)
        assert(problem.G.shape[0] == pi.size)

        # Network quantities
        net.set_var_values(x[:net.num_vars])

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam,nu,mu,pi)
        
