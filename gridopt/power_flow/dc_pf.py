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

class DCPF(PFmethod):
    """
    DC power flow method.
    """

    name = 'DCPF'
    
    _parameters = {'quiet': False,
                   'solver' : 'superlu'}
    
    def __init__(self):

        PFmethod.__init__(self)

        self._parameters = DCPF._parameters.copy()
        self._parameters['solver_parameters'] = {'superlu': {},
                                                 'mumps': {}}

    def create_problem(self,net):

        import pfnet

        # Parameters
        params = self._parameters
        
        # Clear flags
        net.clear_flags()

        # Voltages angles (not slack)
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')

        # Gen active powers (slack)
        net.set_flags('generator',
                      'variable',
                      'slack',
                      'active power')

        # Check
        try:
            assert(net.num_vars == net.num_buses-1+net.get_num_slack_gens())
        except AssertionError:
            raise PFmethodError_BadProblem()

        # Set up problem
        problem = pfnet.Problem(net)
        problem.add_constraint(pfnet.Constraint('DC power balance',net))
        problem.add_constraint(pfnet.Constraint('generator active power participation',net))
        problem.analyze()

        # Return
        return problem
                    
    def solve(self,net):

        from optalg.lin_solver import new_linsolver
        
        # Parameters
        params = self._parameters
        solver_name = params['solver']

        # Copy network
        net = net.get_copy()
        
        # Problem
        t0 = time.time()
        problem = self.create_problem(net)
        problem_time = time.time()-t0
        
        A = problem.A
        b = problem.b
        x = problem.x

        # Solve
        update = True
        t0 = time.time()
        try:
            assert(A.shape[0] == A.shape[1])
            linsolver = new_linsolver(solver_name,'unsymmetric')
            linsolver.analyze(A)
            x = linsolver.factorize_and_solve(A,b)
        except Exception as e:
            update = False
            raise PFmethodError_SolverError(e)
        finally:
            
            # Update network
            if update:
                net.set_var_values(x)
                net.update_properties()
                net.clear_sensitivities()

            # Save results
            self.set_solver_name(solver_name)
            self.set_solver_status('solved' if update else 'error')
            self.set_solver_message('')
            self.set_solver_iterations(1)
            self.set_solver_time(time.time()-t0)
            self.set_solver_primal_variables(x)
            self.set_solver_dual_variables(4*[None])
            self.set_problem(None) # skip for now
            self.set_problem_time(problem_time)
            self.set_network_snapshot(net)
