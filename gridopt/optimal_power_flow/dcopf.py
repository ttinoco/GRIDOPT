#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from scipy.sparse import triu
from method_error import *
from method import OPFmethod
from optalg.opt_solver import OptSolverError,OptSolverIQP,QuadProblem

class DCOPF(OPFmethod):
    """
    DC optimal power flow class.
    """
        
    parameters = {'quiet' : False}
                                    
    def __init__(self):

        OPFmethod.__init__(self)
        self.parameters = DCOPF.parameters.copy()
        self.parameters.update(OptSolverIQP.parameters)

    def create_problem(self,net):
        
        # Parameters
        params = self.parameters
        
        # Clear flags
        net.clear_flags()
        
        # Set flags
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VANG)
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_P_ADJUST,
                      pfnet.GEN_VAR_P)

        # Set up problem
        x = net.get_var_values()
        l = net.get_var_values(pfnet.LOWER_LIMITS)
        u = net.get_var_values(pfnet.UPPER_LIMITS)
        dx = u-l

        constr = pfnet.Constraint(pfnet.CONSTR_TYPE_DCPF,net)
        obj = pfnet.Function(pfnet.FUNC_TYPE_GEN_COST,1.,net)
        
        constr.analyze()
        obj.analyze()
        constr.eval(x)
        obj.eval(x)
        
        H = obj.Hphi + obj.Hphi.T - triu(obj.Hphi)
        g = obj.gphi - H*x

        A = constr.A.copy()
        b = constr.b.copy()

        n = net.num_vars

        try:
            assert(np.all(l < u))
            assert(net.num_vars == (net.num_buses-net.get_num_slack_buses()+
                                    net.get_num_P_adjust_gens()))
            assert(np.abs(obj.phi-(0.5*np.dot(x,H*x)+np.dot(g,x))) < 1e-10)
            assert(H.shape == (n,n))
            assert(A.shape == (net.num_buses,n))
            assert(constr.f.shape == (0,))
            assert(constr.J.shape == (0,n))
        except AssertionError:
            raise OPFmethodError_BadProblem(self)

        # Return
        return QuadProblem(H,g,A,b,l,u)
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters

        # Problem
        problem = self.create_problem(net)
                    
        # Set up solver
        solver = OptSolverIQP()
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError,e:
            raise OPFmethodError_SolverError(self,e)
        finally:
            
            # Get results
            self.results = {'status': solver.get_status(),
                            'error_msg': solver.get_error_msg(),
                            'variables': solver.get_primal_variables(),
                            'iterations': solver.get_iterations()}
            net.update_properties(solver.get_primal_variables())
            self.results.update(net.get_properties())
