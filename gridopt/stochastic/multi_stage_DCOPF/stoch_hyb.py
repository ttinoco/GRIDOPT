#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from method import MS_DCOPF_Method
from problem import MS_DCOPF_Problem
from optalg.stoch_solver import MultiStage_StochHybrid

class MS_DCOPF_SH(MS_DCOPF_Method):
    """
    Stochatic Hybrid Approximation method for multi-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):

        MS_DCOPF_Method.__init__(self)
        self.parameters = MS_DCOPF_SH.parameters.copy()
        self.parameters.update(MS_DCOPF_Problem.parameters)
        self.parameters.update(MultiStage_StochHybrid.parameters)

    def create_problem(self,net,forecast):
        
        problem = MS_DCOPF_Problem(net,forecast)
        problem.set_parameters(self.parameters)
        return problem
        
    def solve(self,net,forecast):
        
        # Local variables
        params = self.parameters
        
        # Parameters
        quiet = params['quiet']

        # Problem
        problem = self.create_problem(net,forecast)
        if not quiet:
            problem.show()

        # Solver
        solver = MultiStage_StochHybrid()
        solver.set_parameters(params)
        
        # Solve
        solver.solve(problem)
        
        
        
