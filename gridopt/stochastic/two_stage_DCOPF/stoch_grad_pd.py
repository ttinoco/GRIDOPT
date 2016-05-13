#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from method import TS_DCOPF_Method
from problem_risk import TS_DCOPF_RiskAverse
from optalg.stoch_solver import PrimalDual_StochGradient

class TS_DCOPF_PDSG(TS_DCOPF_Method):
    """
    Priaml-Dual Stochatic Gradient method for two-stage 
    risk-averse DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):

        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_PDSG.parameters.copy()

    def create_problem(self,net):
        
        return TS_DCOPF_RiskAverse(net)
        
    def solve(self,net):
        
        pass

