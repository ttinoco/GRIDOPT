#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from .method import TS_DCOPF_Method
from .problem_risk import TS_DCOPF_Problem
from .problem_risk import TS_DCOPF_RA_Problem
from optalg.stoch_solver import StochGradientPD

class TS_DCOPF_PDSG(TS_DCOPF_Method):
    """
    Priaml-Dual Stochatic Gradient method for two-stage 
    risk-averse DC OPF problem.
    """
    
    parameters = {'quiet': False}

    name = 'Primal-Dual Stochatic Gradient'
    
    def __init__(self):
        
        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_PDSG.parameters.copy()
        self.parameters.update(TS_DCOPF_Problem.parameters)
        self.parameters.update(TS_DCOPF_RA_Problem.parameters)
        self.parameters.update(StochGradientPD.parameters)

        self.problem = None
        self.results = None

    def create_problem(self,net,parameters):
        
        return TS_DCOPF_RA_Problem(net,parameters)
        
    def solve(self,net):
        
        # Local variables
        params = self.parameters
        
        # Parameters
        quiet = params['quiet']
        
        # Problem
        self.problem = self.create_problem(net,params)
        if not quiet:
            self.problem.show()
            
        # Solver
        solver = StochGradientPD()
        solver.set_parameters(params)
        
        # Solve
        solver.solve(self.problem)

        # Results
        self.results = solver.get_results()
