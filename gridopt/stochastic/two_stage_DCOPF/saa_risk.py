#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import time
import numpy as np
from .method import TS_DCOPF_Method
from .problem import TS_DCOPF_Problem
from .problem_risk import TS_DCOPF_RA_Problem
from scipy.sparse import eye,coo_matrix,bmat

class TS_DCOPF_SAA_Risk(TS_DCOPF_Method):
    """
    SAA method for solving two-stage DCOPF problems with risk constraint.
    """
    
    parameters = {'scenarios': 100,
                  'quiet': False}
    
    def __init__(self):

        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_SAA_Risk.parameters.copy()
        self.parameters.update(TS_DCOPF_Problem.parameters)
        self.parameters.update(TS_DCOPF_RA_Problem.parameters)

        self.problem = None
        self.results = None

    def get_name(self):

        return 'SAA Direct %d' %self.parameters['scenarios']

    def create_problem(self,net,parameters):

        return TS_DCOPF_RA_Problem(net,parameters)

    def solve(self,net):

        # Import
        import cvxpy as cpy

        # Parameters
        params = self.parameters
        num_sce = params['scenarios']
        quiet = params['quiet']
        gamma = params['gamma']

        # Problem
        problemRA = self.create_problem(net,params)
        problem = problemRA.ts_dcopf
        self.problem = problemRA
        if not quiet:
            problemRA.show()

        # Scenarios
        scenarios = [problem.sample_w() for i in range(num_sce)]

        # Constants
        num_p = problem.num_p
        num_w = problem.num_w
        num_r = problem.num_r
        num_bus = problem.num_bus
        num_br = problem.num_br

        H0 = problem.H0
        g0 = problem.g0
        H1 = problem.H1
        g1 = problem.g1
        G = problem.G
        A = problem.A
        b = problem.b
        R = problem.R
        J = problem.J
        p_min = problem.p_min
        p_max = problem.p_max
        t_min = params['t_min']*problemRA.Qref
        t_max = params['t_max']*problemRA.Qref
        z_max = problem.z_max
        z_min = problem.z_min

        # Variables
        p = cpy.Variable(num_p)
        t = cpy.Variable(1)
        q_list = [cpy.Variable(num_p) for i in range(num_sce)]
        w_list = [cpy.Variable(num_w) for i in range(num_sce)]
        s_list = [cpy.Variable(num_r) for i in range(num_sce)]

        # Objective
        obj = 0.5*cpy.quad_form(p,cpy.Constant(H0)) + p.T*cpy.Constant(g0)
        for i in range(num_sce):
            obj = obj + (1./float(num_sce))*(0.5*cpy.quad_form(q_list[i],cpy.Constant(H1)) + q_list[i].T*cpy.Constant(g1))

        # Constraints
        constr = [p_min <= p, p <= p_max, t_min <= t, t <= t_max]
        expr = (1.-gamma)*t
        for i in range(num_sce):
            expr = expr + (1./float(num_sce))*cpy.max_elemwise(0.5*cpy.quad_form(q_list[i],cpy.Constant(H1)) + q_list[i].T*cpy.Constant(g1) - 
                                                               problemRA.Qmax - t,0)
        constr.append(expr <= 0)
        for i in range(num_sce):
            constr.append(p_min <= p + q_list[i])
            constr.append(p + q_list[i] <= p_max)
            constr.append(0 <= s_list[i])
            constr.append(s_list[i] <= scenarios[i])
            constr.append(z_min <= cpy.Constant(J)*w_list[i])
            constr.append(cpy.Constant(J)*w_list[i] <= z_max)
            constr.append(cpy.Constant(G)*(p+q_list[i]) + 
                          cpy.Constant(R)*s_list[i] -
                          cpy.Constant(A)*w_list[i] == b)
            
        prob = cpy.Problem(cpy.Minimize(obj),constr)
        prob.solve(solver=cpy.ECOS,verbose=True)
        self.results = [(0,0.,np.hstack((np.array(p.value).flatten(),np.array(t.value).flatten())),np.nan)]
        
        
        
