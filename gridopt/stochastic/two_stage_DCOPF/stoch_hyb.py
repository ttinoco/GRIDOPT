#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from problem import TS_DCOPF
from method import TS_DCOPF_Method
from optalg.stoch_solver import StochHybrid

class TS_DCOPF_SH(TS_DCOPF_Method):
    """
    Stochatic Hybrid Approximation method for two-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):

        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_SH.parameters.copy()

    def create_problem(self,net):
        
        return TS_DCOPF(net)
        
    def solve(self,net):
        
        pass

