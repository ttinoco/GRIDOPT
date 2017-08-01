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
from numpy.linalg import norm

class DCOPF(PFmethod):
    """
    DC optimal power flow method.
    """

    name = 'DCOPF'
        
    parameters = {'thermal_limits': False,
                  'renewable_curtailment': False,
                  'optsolver' : 'iqp'}

    parameters_iqp = {}

    parameters_augl = {}

    parameters_ipopt = {}
    
    def __init__(self):

        from optalg.opt_solver import OptSolverIQP, OptSolverAugL, OptSolverIpopt

        # Parent init
        PFmethod.__init__(self)
        
        # Optsolver parameters
        iqp_params = OptSolverIQP.parameters.copy()
        iqp_params.update(self.parameters_iqp)     # overwrite defaults
        augl_params = OptSolverAugL.parameters.copy()
        augl_params.update(self.parameters_augl)   # overwrite defaults
        ipopt_params = OptSolverIpopt.parameters.copy()
        ipopt_params.update(self.parameters_ipopt) # overwrite defaults
        self.parameters.update(DCOPF.parameters)
        self.parameters['optsolver_params'] = {'iqp': iqp_params,
                                               'augl': augl_params,
                                               'ipopt': ipopt_params}

    def create_problem(self,net):

        import pfnet
        
        # Parameters
        params = self.parameters
        thermal_limits = params['thermal_limits']
        
        # Clear flags
        net.clear_flags()
        
        # Set flags
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')
        net.set_flags('generator',
                      ['variable','bounded'],
                      ['adjustable active power','not on outage'],
                      'active power')
        net.set_flags('load',
                      ['variable','bounded'],
                      'adjustable active power',
                      'active power')
        if params['renewable_curtailment']:
            net.set_flags('variable generator',
                          ['variable','bounded'],
                          'any',
                          'active power')

        try:
            num_gvar =  len([g for g in net.generators if 
                             (not g.is_on_outage()) and g.is_P_adjustable()])
            num_cur = net.num_var_generators if params['renewable_curtailment'] else 0
            assert(net.num_bounded == (num_gvar+net.get_num_P_adjust_loads()+num_cur)*net.num_periods)
            assert(net.num_vars == (net.num_buses-net.get_num_slack_buses()+
                                    num_gvar+net.get_num_P_adjust_loads()+
                                    num_cur)*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
            
        # Set up problem
        problem = pfnet.Problem(net)
        problem.add_constraint(pfnet.Constraint('variable bounds',net))
        problem.add_constraint(pfnet.Constraint('DC power balance',net))
        if thermal_limits:
            problem.add_constraint(pfnet.Constraint('DC branch flow limits',net))
        problem.add_function(pfnet.Function('generation cost',1.,net))
        problem.add_function(pfnet.Function('consumption utility',-1.,net))
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net):

        from optalg.opt_solver import OptSolverError
        from optalg.opt_solver import OptSolverIQP, OptSolverAugL, OptSolverIpopt
        
        # Parameters
        params = self.parameters
        optsolver_name = params['optsolver']
        optsolver_params = params['optsolver_params']

        # Opt solver
        if optsolver_name == 'iqp':
            optsolver = OptSolverIQP()
        elif optsolver_name == 'augl':
            optsolver = OptSolverAugL()
        elif optsolver_name == 'ipopt':
            optsolver = OptSolverIpopt()
        else:
            raise PFmethodError_BadOptSolver(self)
        optsolver.set_parameters(optsolver_params[optsolver_name])

        # Problem
        problem = self.create_problem(net)
        
        # Set up solver
        solver = OptSolverIQP()
        solver.set_parameters(params)
        
        # Solve
        try:
            optsolver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:

            # Update net properties
            net.update_properties(optsolver.get_primal_variables()[:net.num_vars])
            
            # Get results
            self.set_status(optsolver.get_status())
            self.set_error_msg(optsolver.get_error_msg())
            self.set_iterations(optsolver.get_iterations())
            self.set_primal_variables(optsolver.get_primal_variables())
            self.set_dual_variables(optsolver.get_dual_variables())
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
        assert(net.num_vars+problem.num_extra_vars == x.size)
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
