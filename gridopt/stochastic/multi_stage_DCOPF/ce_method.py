#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import time
import numpy as np
from problem import MS_DCOPF
from method import MS_DCOPF_Method
from scipy.sparse import eye,coo_matrix,bmat
from optalg.opt_solver import OptSolverIQP, QuadProblem

class MS_DCOPF_ce(MS_DCOPF_Method):
    """
    Certainty-Equivalent metohd for multi-stage DC OPF problem.
    """
    
    parameters = {'quiet': False}
    
    def __init__(self):

        MS_DCOPF_Method.__init__(self)
        self.parameters = MC_DCOPF_ce.parameters.copy()

    def create_problem(self,net,forecast):
        
        return None
        
    def solve(self,net,forecast):
        
        pass

