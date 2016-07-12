#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from method import MS_DCOPF_Method
from problem import MS_DCOPF_Problem
from optalg.stoch_solver import StochDualDynProg,StochProblemMS_Tree

class MS_DCOPF_SDDP(MS_DCOPF_Method):
    """
    Stochastic dual dynamic programming method for multi-stage DC OPF problem.
    """
    
    parameters = {'branching_factor': 2,
                  'quiet': False}
    
    def __init__(self):
        
        MS_DCOPF_Method.__init__(self)
        self.parameters = MS_DCOPF_SDDP.parameters.copy()
        self.parameters.update(MS_DCOPF_Problem.parameters)
        self.parameters.update(StochDualDynProg.parameters)

    def create_problem(self,net,forecast,parameters):
        
        return MS_DCOPF_Problem(net,forecast,parameters)
        
    def solve(self,net,forecast):
        
        # Local variables
        params = self.parameters
        
        # Parameters
        quiet = params['quiet']
        bfactor = params['branching_factor']

        # Problem
        self.problem = self.create_problem(net,forecast,params)

        # Scenario tree
        self.tree = StochProblemMS_Tree(self.problem,bfactor)

        # Show
        if not quiet:
            self.tree.show()
            self.problem.show(scenario_tree=self.tree)
   
        # Solver
        solver = StochDualDynProg()
        solver.set_parameters(params)
        
        # Solve
        solver.solve(self.problem,self.tree)
        
        # Return policy
        return solver.get_policy()
