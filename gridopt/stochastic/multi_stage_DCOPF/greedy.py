#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import time
import numpy as np
from types import MethodType
from method import MS_DCOPF_Method
from problem import MS_DCOPF_Problem
from scipy.sparse import eye,coo_matrix,bmat
from optalg.stoch_solver import StochObjMS_Policy
from optalg.opt_solver import OptSolverIQP, QuadProblem

class MS_DCOPF_GR(MS_DCOPF_Method):
    """
    Greedy method for multi-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):
        
        MS_DCOPF_Method.__init__(self)
        self.parameters = MS_DCOPF_GR.parameters.copy()

    def create_problem(self,net,forecast):
        
        return MS_DCOPF_Problem(net,forecast)
        
    def solve(self,net,forecast):
        
        # Local variables
        params = self.parameters

        # Parameters
        quiet = params['quiet']

        # Problem
        self.problem = self.create_problem(net,forecast)
        if not quiet:
            self.problem.show()
 
        # Construct policy
        def apply(cls,t,x_prev,Wt):
            
            assert(0 <= t < cls.problem.T)
            assert(len(Wt) == t+1)
            
            x_list,Q_list,gQ_list = cls.problem.eval_stage_approx(t,[Wt[-1]],x_prev,quiet=True,tf=t)
            assert(len(x_list) == 1)
            x = x_list[0]
            p,q,w,s,y,z = cls.problem.separate_x(x)
            p_prev = x_prev[:cls.problem.num_p]
            
            # Check dims
            assert(p.shape,(cls.problem.num_p,))
            assert(p_prev.shape,(cls.problem.num_p,))
            
            # Check feasibility
            if not cls.problem.is_point_feasible(t,p,p_prev,q,w,s,z,Wt[-1]):
                raise ValueError('infeasible point')
            
            # Return
            return cls.problem.construct_x(p=p,q=q,w=w,s=s,y=p-p_prev,z=z)
            
        policy = StochObjMS_Policy(self.problem,data=None,name='Greedy')
        policy.apply = MethodType(apply,policy)
        
        # Return
        return policy
