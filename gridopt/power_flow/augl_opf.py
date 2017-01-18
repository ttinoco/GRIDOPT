#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import pfnet
import numpy as np
from .method_error import *
from .method import PFmethod
from optalg.opt_solver import OptSolverError,OptCallback,OptTermination,OptSolverAugL

class AugLOPF(PFmethod):
    """
    Augmented Lagrangian-based optimal power flow method.
    """
 
    name = 'AugLOPF'

    parameters = {'weight_cost':1e-2,   # for generation cost
                  'weight_limit':1e2,   # for soft limits
                  'weight_reg':1e-5,    # for regularization
                  'vmin_thresh':0.1}    # threshold for vmin
                   
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = AugLOPF.parameters.copy()
        self.parameters.update(OptSolverAugL.parameters)

    def create_problem(self,net):

        # Parameters
        params = self.parameters
        wc = params['weight_cost']
        wl = params['weight_limit']
        wr = params['weight_reg']
        
        # Clear flags
        net.clear_flags()
        
        # Voltage magnitudes
        net.set_flags('bus',
                      'variable',
                      'any',
                      'voltage magnitude')
        
        # Voltage angles
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')

        # Generator active power
        net.set_flags('generator',
                      ['variable','bounded'],
                      'not on outage',
                      'active power')

        # Generator reactive power
        net.set_flags('generator',
                      ['variable','bounded'],
                      'regulator',
                      'reactive power')

        try:
            assert(net.num_vars == (2*net.num_buses-net.get_num_slack_buses() +
                                    net.get_num_gens_not_on_outage() + 
                                    net.get_num_reg_gens())*net.num_periods)
            assert(net.num_bounded == (net.get_num_gens_not_on_outage() + 
                                       net.get_num_reg_gens())*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
                                    
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint('AC power balance')
        problem.add_constraint('variable nonlinear bounds') 
        problem.add_function('generation cost',wc)
        problem.add_function('soft voltage magnitude limits',wl)
        problem.add_function('voltage angle regularization',wr)
        problem.add_function('generator powers regularization',wr)
        problem.analyze()
        
        # Return
        return problem

    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.network
            if header:
                print('{0:^5}'.format('vmax'), end=' ')
                print('{0:^5}'.format('vmin'), end=' ')
                print('{0:^6}'.format('bvvio'), end=' ')
                print('{0:^6}'.format('gQvio'), end=' ')
                print('{0:^6}'.format('gPvio'))
            else:
                print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.bus_v_vio)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.gen_Q_vio)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.gen_P_vio)))
        return info_printer
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        vmin_thresh = params['vmin_thresh']

        # Problem
        problem = self.create_problem(net)

        # Callbacks

        # Termination
        def t1(s):
            if np.min(s.problem.network.bus_v_min) < vmin_thresh:
                return True
            else:
                return False
        
        # Info printer
        info_printer = self.get_info_printer()
        
        # Set up solver
        solver = OptSolverAugL()
        solver.set_parameters(params)
        solver.add_termination(OptTermination(t1,'low voltage'))
        solver.set_info_printer(info_printer)
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:
            
            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(solver.get_primal_variables())
            self.set_dual_variables(solver.get_dual_variables())
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)

    def update_network(self,net):
        
        # Get data
        problem = self.results['problem']
        x = self.results['primal_variables']
        lam,nu,mu,pi = self.results['dual_variables']
       
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == x.shape)
        assert(net.num_vars == x.size)
        assert(problem.A.shape[0] == lam.size)
        assert(problem.f.shape[0] == nu.size)
        assert(mu is None or not mu.size)
        assert(pi is None or not pi.size)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam,nu,mu,pi)
        
