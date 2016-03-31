#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from method_error import *
from method import PFmethod
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
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VANG)

        # Gen active powers (slack)
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_SLACK,
                      pfnet.GEN_VAR_P)

        # Check
        try:
            assert(net.num_vars == net.num_buses-1+net.get_num_slack_gens())
        except AssertionError:
            raise PFmethodError_BadProblem(self)

        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_DCPF)
        problem.add_constraint(pfnet.CONSTR_TYPE_PAR_GEN_P)
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
        except Exception,e:
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
        assert(lam is None)
        assert(nu is None)
        assert(mu is None)
        assert(pi is None)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
