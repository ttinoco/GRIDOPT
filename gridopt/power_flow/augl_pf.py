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

class AugLPF(PFmethod):
    """
    Augmented Lagrangian-based power flow method.
    """

    name = 'AugLPF'
        
    parameters = {'weight_vmag': 1e0,  # for reg voltage magnitude penalty
                  'weight_vang': 1e0,  # for angle difference penalty
                  'weight_pq': 1e-3,   # for gen powers
                  'weight_t': 1e-3,    # for tap ratios
                  'weight_b': 1e-3,    # for shunt susceptances
                  'lock_taps': True,   # flag for locking transformer tap ratios
                  'lock_shunts': True, # flag for locking swtiched shunts
                  'feastol' : 1e-4,    # see AugL
                  'optol' : 1e-4,      # see AugL
                  'kappa' : 1e-4,      # see AugL
                  'vmin_thresh': 0.1}  # threshold for vmin
                  
    def __init__(self):

        from optalg.opt_solver import OptSolverAugL

        PFmethod.__init__(self)
        parameters = OptSolverAugL.parameters.copy()
        parameters.update(AugLPF.parameters)
        self.parameters = parameters

    def create_problem(self,net):

        import pfnet

        # Parameters
        params = self.parameters
        wm = params['weight_vang']
        wa = params['weight_vmag']
        wp = params['weight_pq']
        wt = params['weight_t']
        wb = params['weight_b']
        lock_taps = params['lock_taps']
        lock_shunts = params['lock_shunts']
        
        # Clear flags
        net.clear_flags()

        # Adjust gens
        net.adjust_generators()
        
        # Set up variables
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      ['voltage magnitude','voltage angle'])
        net.set_flags('bus',
                      'variable',
                      ['not slack','regulated by generator'],
                      'voltage magnitude deviation')
        net.set_flags('generator',
                      'variable',
                      'slack',
                      'active power')
        net.set_flags('generator',
                      'variable',
                      'regulator',
                      'reactive power')
        
        # Tap ratios
        if not lock_taps:
            net.set_flags('bus',
                          'variable',
                          'regulated by transformer',
                          'voltage magnitude violation')
            net.set_flags('branch', 
                          'variable',
                          'tap changer - v',
                          ['tap ratio','tap ratio deviation'])

        # Shunt voltage control
        if not lock_shunts:
            net.set_flags('bus',
                          'variable',
                          'regulated by shunt',
                          'voltage magnitude violation')
            net.set_flags('shunt', 
                          'variable',
                          'switching - v',
                          ['susceptance','susceptance deviation'])

        try:
            num_vars = (2*(net.num_buses-net.get_num_slack_buses()) +
                        2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                        net.get_num_slack_gens() +
                        net.get_num_reg_gens())*net.num_periods
            if not lock_taps:
                num_vars += (2*net.get_num_buses_reg_by_tran() +
                             3*net.get_num_tap_changers_v())*net.num_periods
            if not lock_shunts:
                num_vars += (2*net.get_num_buses_reg_by_shunt() +
                             3*net.get_num_switched_shunts())*net.num_periods
            if not lock_shunts and not lock_taps:
                num_vars -= 2*len([b for b in net.buses 
                                   if (b.is_regulated_by_tran() and
                                       b.is_regulated_by_shunt())])*net.num_periods
            assert(net.num_vars == num_vars)
        except AssertionError:
            raise PFmethodError_BadProblem(self)  
            
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.Constraint('AC power balance',net))
        problem.add_constraint(pfnet.Constraint('generator active power participation',net))
        problem.add_constraint(pfnet.Constraint('generator reactive power participation',net))
        problem.add_constraint(pfnet.Constraint('voltage regulation by generators',net))
        if not lock_taps:
            problem.add_constraint(pfnet.Constraint('voltage regulation by transformers',net))
        if not lock_shunts:
            problem.add_constraint(pfnet.Constraint('voltage regulation by shunts',net))
        problem.add_function(pfnet.Function('voltage magnitude regularization',wm/max([net.num_buses,1.]),net))
        problem.add_function(pfnet.Function('voltage angle regularization',wa/max([net.num_buses,1.]),net))
        problem.add_function(pfnet.Function('generator powers regularization',wp/max([net.num_generators,1.]),net))
        if not lock_taps:
            problem.add_function(pfnet.Function('tap ratio regularization',wt/max([net.get_num_tap_changers_v(),1.]),net))
        if not lock_shunts:
            problem.add_function(pfnet.Function('susceptance regularization',wb/max([net.get_num_switched_shunts(),1.]),net))
        problem.analyze()
        
        # Return
        return problem

    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.network
            if header:
                print('{0:^5}'.format('vmax'), end=' ')
                print('{0:^5}'.format('vmin'), end=' ')
                print('{0:^8}'.format('gvdev'), end=' ')
                print('{0:^8}'.format('gQvio'), end=' ')
                print('{0:^8}'.format('tvvio'), end=' ')
                print('{0:^8}'.format('svvio'))
            else:
                print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.gen_v_dev)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.gen_Q_vio)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.tran_v_vio)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.shunt_v_vio)))
        return info_printer
            
    def solve(self,net):

        from optalg.opt_solver import OptSolverError, OptCallback, OptTermination, OptSolverAugL
        
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
           
            # Restore net properties
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
        assert(problem.x.size == mu.size)
        assert(problem.x.size == pi.size)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam,nu,np.zeros(0),np.zeros(0))
