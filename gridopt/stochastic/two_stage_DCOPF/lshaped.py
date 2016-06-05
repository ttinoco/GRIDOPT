#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import time
import numpy as np
from .method import TS_DCOPF_Method
from .problem import TS_DCOPF_Problem
from scipy.sparse import eye,coo_matrix,bmat
from optalg.opt_solver import OptSolverIQP, QuadProblem

class TS_DCOPF_LS(TS_DCOPF_Method):
    """
    L-shaped method for solving two-stage DCOPF problems.
    """
    
    parameters = {'scenarios':100,
                  'maxiters':30,
                  'samples':300,
                  't_max':1e6,
                  'z_max':1e6}
    
    def __init__(self):

        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_LS.parameters.copy()
    
    def create_problem(self,net):

        return TS_DCOPF_Problem(net)

    def solve(self,net):

        # Parameters
        params = self.parameters
        num_sce = params['scenarios']
        maxiters = params['maxiters']
        samples = params['samples']
        t_max = params['t_max']
        z_max = params['z_max']

        # Problem
        problem = self.create_problem(net)

        # Initial point
        p,results = problem.solve_approx(quiet=True)

        # Scenarios
        scenarios = [problem.sample_w() for i in range(num_sce)]

        # Matvecs
        H0 = problem.H0
        g0 = problem.g0
        p_min = problem.p_min
        p_max = problem.p_max
        A1 = np.zeros((0,p.size+1))
        b = np.zeros(0)

        # Init eval
        Q,gQ = [sum(l)/float(num_sce) for l in zip(*[problem.eval_Q(p,r) for r in scenarios])]

        # Iteartions
        t = 0
        num_z = 0
        t0 = time.time()
        solved = False
        for k in range(maxiters):

            # Show info
            t1 = time.time()
            F = 0.5*np.dot(p,H0*p)+np.dot(g0,p) + Q
            EF,EgF = problem.eval_EF(p,samples=samples)
            print('%d,%.2f,%.2e,%.2e,%.5e,%.5e' %(k,t1-t0,Q,t,F,EF))
            t0 += time.time()-t1
            
            # Solved
            if solved:
                print('solved')
                break
            
            # Add cut
            a = np.hstack((-gQ,1.))
            A1 = np.vstack((A1,a))
            b = np.hstack((b,Q-np.dot(gQ,p)))
            num_z += 1
            
            # Prepare problem
            A1sp = coo_matrix(A1)
            A = bmat([[A1sp,-eye(num_z)]],format='coo')
            
            H = bmat([[H0,None,None],
                      [None,coo_matrix((1,1)),None],
                      [None,None,coo_matrix((num_z,num_z))]],format='coo')

            g = np.hstack((g0,1,np.zeros(num_z)))
            
            u = np.hstack((p_max,t_max,z_max*np.ones(num_z)))
            l = np.hstack((p_min,-t_max,np.zeros(num_z)))
            qp_prob = QuadProblem(H,g,A,b,l,u)

            # Solve problem
            qp_sol = OptSolverIQP()
            qp_sol.set_parameters({'quiet':True})
            qp_sol.solve(qp_prob)

            # Get results
            x = qp_sol.get_primal_variables()
            p = x[:p.size]
            t = x[p.size]
            z = x[p.size+1:]
            assert(np.all(z > 0))
            
            # Eval
            Q,gQ = [sum(l)/float(num_sce) for l in zip(*[problem.eval_Q(p,r) for r in scenarios])]
            
            # Check solved
            if Q <= t:
                solved = True
            
            
            
            
            
            

        
        
