#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from method import TwoStageDCOPF_Method
from problem import TwoStageDCOPF_Problem

class LShaped(TwoStageDCOPF_Method):
    """
    L-shaped method for solving two-stage DCOPF problems.
    """

    parameters = {'scenarios':50}

    def __init__(self):

        TwoStageDCOPF_Method.__init__(self)
        self.parameters = LShaped.parameters.copy()
    
    def create_problem(self,net):

        return TwoStageDCOPF_Problem(net)

    def solve(self,net):

        # Problem
        problem = self.create_problem(net)

        # Initial point
        x_ce = problem.solve_certainty_equivalent()

        
        

        
        
