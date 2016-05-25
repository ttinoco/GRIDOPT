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

class MS_DCOPF_TSCE(MS_DCOPF_Method):
    """
    Certainty-Equivalent metohd for multi-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):
        
        MS_DCOPF_Method.__init__(self)
        self.parameters = MS_DCOPF_TSCE.parameters.copy()
        self.parameters.update(MS_DCOPF_Problem.parameters)
        
    def create_problem(self,net,forecast,parameters):
        
        return MS_DCOPF_Problem(net,forecast,parameters) 
        
    def solve(self,net,forecast):
        
        # Local variables
        params = self.parameters

        # Parameters
        quiet = params['quiet']

        # Problem
        self.problem = self.create_problem(net,forecast,params)
        if not quiet:
            self.problem.show()

        # Prediction
        Er_list = self.problem.predict_W(self.problem.get_num_stages()-1)
        assert(len(Er_list) == self.problem.T)

        # Solve certainty equivalent
        x_list,Q_list,gQ_list,results = self.problem.eval_stage_approx(0,Er_list,self.problem.x_prev,quiet=quiet)
        
        # Slow generator powers
        p_list = [x[:self.problem.num_p] for x in x_list]
        assert(len(p_list) == self.problem.T)
 
        # Construct policy
        def apply(cls,t,x_prev,Wt):
            
            assert(0 <= t < cls.problem.T)
            assert(len(Wt) == t+1)
            
            p = cls.data[t]
            p_prev = x_prev[:cls.problem.num_p]
            
            q,w,s,z = cls.problem.eval_stage_adjust(t,Wt[-1],p,quiet=True)

            x = cls.problem.construct_x(p=p,q=q,w=w,s=s,y=p-p_prev,z=z)

            # Check feasibility
            if not cls.problem.is_point_feasible(t,x,x_prev,Wt[-1]):
                raise ValueError('point not feasible')
            
            # Return
            return x
            
        policy = StochObjMS_Policy(self.problem,data=p_list,name='Two-Stage Certainty-Equivalent')
        policy.apply = MethodType(apply,policy)
        
        # Return
        return policy
