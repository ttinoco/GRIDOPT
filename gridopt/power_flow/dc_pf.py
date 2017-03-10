#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from .method_error import *
from .method import PFmethod
from optalg.lin_solver import new_linsolver

class DCPF(PFmethod):
    """
    DC power flow method.
    """

    name = 'DCPF'
    
    parameters = {'quiet': False}
    
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = DCPF.parameters.copy()

    def create_problem(self,net):

        # Parameters
        params = self.parameters
        
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
            raise PFmethodError_BadProblem(self)

        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint('DC power balance')
        problem.add_constraint('generator active power participation')
        problem.analyze()

        # Return
        return problem
                    
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        
        # Problem
        problem = self.create_problem(net)
        A = problem.A
        b = problem.b
        x = problem.x

        # Solve
        try:
            assert(A.shape[0] == A.shape[1])
            linsolver = new_linsolver('default','unsymmetric')
            linsolver.analyze(A)
            x = linsolver.factorize_and_solve(A,b)
            net.update_properties(x)
            self.set_status('solved')
        except Exception as e:
            raise PFmethodError_SolverError(self,e)
        finally:
            
            # Update net properties
            net.update_properties(x)

            # Get results
            self.set_iterations(1)
            self.set_primal_variables(x)
            self.set_dual_variables(4*[None])
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
        assert(lam is None or not lam.size)
        assert(nu is None or not nu.size)
        assert(mu is None or not mu.size)
        assert(pi is None or not pi.size)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
