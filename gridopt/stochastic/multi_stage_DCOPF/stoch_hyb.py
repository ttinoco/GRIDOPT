#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from .method import MS_DCOPF_Method
from .problem import MS_DCOPF_Problem
from optalg.stoch_solver import StochHybridMS

class MS_DCOPF_SH(MS_DCOPF_Method):
    """
    Stochatic Hybrid Approximation method for multi-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):
        
        MS_DCOPF_Method.__init__(self)
        self.parameters = MS_DCOPF_SH.parameters.copy()
        self.parameters.update(MS_DCOPF_Problem.parameters)
        self.parameters.update(StochHybridMS.parameters)

    def create_problem(self,net,forecast,parameters):
        
        return MS_DCOPF_Problem(net,forecast,parameters)
        
    def solve(self,net,forecast):
        
        # Local variables
        params = self.parameters
        
        # Parameters
        quiet = params['quiet']

        # Problem
        problem = self.create_problem(net,forecast,params)
        if not quiet:
            problem.show()
            
        # Solver
        solver = StochHybridMS()
        solver.set_parameters(params)
        
        # Solve
        solver.solve(problem)
        
        # Return policy
        return solver.get_policy()
        
