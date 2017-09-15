#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import time
import numpy as np
from .method_error import *
from .method import PFmethod
from numpy.linalg import norm

class DCOPF(PFmethod):
    """
    DC optimal power flow method.
    """

    name = 'DCOPF'
        
    _parameters = {'thermal_limits': False,
                   'renewable_curtailment': False,
                   'optsolver' : 'iqp'}

    _parameters_iqp = {}

    _parameters_augl = {}

    _parameters_ipopt = {}
    
    def __init__(self):

        from optalg.opt_solver import OptSolverIQP, OptSolverAugL, OptSolverIpopt

        # Parent init
        PFmethod.__init__(self)
        
        # Optsolver parameters
        iqp_params = OptSolverIQP.parameters.copy()
        iqp_params.update(self._parameters_iqp)     # overwrite defaults

        augl_params = OptSolverAugL.parameters.copy()
        augl_params.update(self._parameters_augl)   # overwrite defaults

        ipopt_params = OptSolverIpopt.parameters.copy()
        ipopt_params.update(self._parameters_ipopt) # overwrite defaults

        self._parameters.update(DCOPF._parameters)
        self._parameters['optsolver_parameters'] = {'iqp': iqp_params,
                                                    'augl': augl_params,
                                                    'ipopt': ipopt_params}

    def create_problem(self,net):

        import pfnet
        
        # Parameters
        params = self._parameters
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
            raise PFmethodError_BadProblem()
            
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
        params = self._parameters
        optsolver_name = params['optsolver']
        optsolver_params = params['optsolver_parameters']

        # Opt solver
        if optsolver_name == 'iqp':
            optsolver = OptSolverIQP()
        elif optsolver_name == 'augl':
            optsolver = OptSolverAugL()
        elif optsolver_name == 'ipopt':
            optsolver = OptSolverIpopt()
        else:
            raise PFmethodError_BadOptSolver()
        optsolver.set_parameters(optsolver_params[optsolver_name])

        # Copy network
        net = net.get_copy()

        # Problem
        t0 = time.time()
        problem = self.create_problem(net)
        problem_time = time.time()-t0
        
        # Set up solver
        solver = OptSolverIQP()
        solver.set_parameters(params)
        
        # Solve
        t0 = time.time()
        try:
            optsolver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(e)
        finally:

            # Update network
            net.set_var_values(optsolver.get_primal_variables()[:net.num_vars])
            net.update_properties()
            net.clear_sensitivities()
            problem.store_sensitivities(*optsolver.get_dual_variables())

            # Save results
            self.set_solver_status(optsolver.get_status())
            self.set_solver_message(optsolver.get_error_msg())
            self.set_solver_iterations(optsolver.get_iterations())
            self.set_solver_time(time.time()-t0)
            self.set_solver_primal_variables(optsolver.get_primal_variables())
            self.set_solver_dual_variables(optsolver.get_dual_variables())
            self.set_problem(None) # skip for now
            self.set_problem_time(problem_time)
            self.set_network_snapshot(net)
