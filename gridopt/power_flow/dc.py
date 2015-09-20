#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from method_error import *
from method import PFmethod
from optalg.lin_solver import LinSolverMUMPS

class DCPF(PFmethod):
    """
    DC power flow method.
    """
    
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
            linsolver = LinSolverMUMPS('unsymmetric')
            linsolver.analyze(A)
            x = linsolver.factorize_and_solve(A,b)
            net.update_properties(x)
            self.set_status('solved')
        except Exception,e:
            raise PFmethodError_SolverError(self,e)
        finally:

            # Get results
            self.results = {'status': self.results['status'],
                            'error_msg': self.results['error_msg'],
                            'variables': x,
                            'iterations': 0}
            self.results.update(net.get_properties())
