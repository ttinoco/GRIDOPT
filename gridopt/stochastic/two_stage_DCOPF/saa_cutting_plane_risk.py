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
from optalg.opt_solver import OptSolverIQP, QuadProblem

class TS_DCOPF_SAA_CP_Risk(TS_DCOPF_Method):
    """
    SAA cutting-plane method for solving two-stage DCOPF problems with risk constraint.
    """
    
    parameters = {'scenarios': 100,
                  'num_procs': 1,
                  'maxiters': 100,
                  'maxtime': 600,
                  'period': 60,
                  'z_inf': 1e6,
                  'y_inf': 1e6,
                  'quiet': False}
    
    def __init__(self):

        TS_DCOPF_Method.__init__(self)
        self.parameters = TS_DCOPF_SAA_CP_Risk.parameters.copy()
        self.parameters.update(TS_DCOPF_Problem.parameters)
        self.parameters.update(TS_DCOPF_RA_Problem.parameters)

        self.problem = None
        self.results = None

    def get_name(self):

        return 'SAA Cutting-Plane %d' %self.parameters['scenarios']

    def create_problem(self,net,parameters):

        return TS_DCOPF_RA_Problem(net,parameters)

    def solve(self,net):

        # Imports
        from multiprocess import Pool

        # Parameters
        params = self.parameters
        num_sce = params['scenarios']
        num_procs = params['num_procs']
        maxiters = params['maxiters']
        maxtime = params['maxtime']
        period = params['period']
        z_inf = params['z_inf']
        y_inf = params['y_inf']
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

        # Pool
        pool = Pool(num_procs)

        # Constants
        num_p = problem.num_p
        H0 = problem.H0
        g0 = problem.g0
        p_min = problem.p_min
        p_max = problem.p_max
        t_min = params['t_min']*problemRA.Qref
        t_max = params['t_max']*problemRA.Qref
        op = np.zeros(num_p)        

        # For base eq
        A0 = bmat([[op,1.-gamma,0,1.,-1.]]) # p t z1 z2 z3
        b0 = np.zeros(1)

        # For obj cuts eq
        A1 = np.zeros((0,num_p+4))
        b1 = np.zeros(0)

        # For constr cuts eq
        A2 = np.zeros((0,num_p+4))
        b2 = np.zeros(0)

        # Header
        if not quiet:
            print('\nSAA Cutting-Plane Risk')
            print('-----------------------')
            print('{0:^8s}'.format('iter'), end=' ')
            print('{0:^10s}'.format('time(s)'), end=' ')
            print('{0:^12s}'.format('FL'), end= ' ')
            print('{0:^12s}'.format('GL'), end= ' ')
            print('{0:^12s}'.format('saved'))

        # Init
        k = 0
        t1 = 0
        num_y = 0
        t0 = time.time()
        self.results = []

        # Loop
        while True:

            # Construct master
            H = bmat([[H0,coo_matrix((num_p,4+2*num_y))],
                      [coo_matrix((4+2*num_y,num_p)),None]],format='coo')
            g = np.hstack((g0,0.,1.,np.zeros(2+2*num_y)))
            if num_y > 0:
                A = bmat([[A0,None,None],
                          [A1,-eye(num_y),None],
                          [A2,None,-eye(num_y)]],format='coo')
            else:
                A = bmat([[A0]],format='coo')
            b = np.hstack((b0,b1,b2))
            l = np.hstack((p_min,  # p
                           t_min,  # t
                           -z_inf, # z1
                           -z_inf, # z2
                           -z_inf, # z3
                           -y_inf*np.ones(num_y),  # y1
                           -y_inf*np.ones(num_y))) # y2
            u = np.hstack((p_max,  # p
                           t_max,  # t
                           z_inf,  # z1
                           z_inf,  # z2
                           0,      # z3
                           np.zeros(num_y),  # y1
                           np.zeros(num_y))) # y2
            qp = QuadProblem(H,g,A,b,l,u)

            # Solve problem
            solver = OptSolverIQP()
            solver.set_parameters({'quiet': True,
                                   'tol': self.parameters['tol']})
            solver.solve(qp)
            assert(solver.get_status() == 'solved')
            x_full = solver.get_primal_variables()
            p = x_full[:num_p]
            t = x_full[num_p]
            z1 = x_full[num_p+1]
            z2 = x_full[num_p+2]
            self.x = x_full[:num_p+1]

            # Save
            if time.time()-t0 > t1:
                self.results.append((k,time.time()-t0,self.x,np.nan))
                t1 += period

            # Iters
            if k >= maxiters:
                break
                
            # Maxtime
            if time.time()-t0 >= maxtime:
                break
            
            # Obj saa
            Q_list,gQ_list = zip(*pool.map(lambda w: problem.eval_Q(p,w),scenarios))
            Q = sum(Q_list)/float(num_sce)
            gQ = sum(gQ_list)/float(num_sce)

            # Constr saa
            S_list = []
            gS_list = []
            for i in range(num_sce):
                qq = Q_list[i]-problemRA.Qmax-t
                S_list.append(np.maximum(qq,0.))
                if qq >= 0.:
                    gS_list.append(np.hstack((gQ_list[i],-1.)))
                else:
                    gS_list.append(np.hstack((op,0.)))
            S = sum(S_list)/float(num_sce)
            gS = sum(gS_list)/float(num_sce)

            # Obj cut
            a1 = np.hstack((gQ,0.,-1.,0.,0.))
            A1 = np.vstack((A1,a1))
            b1 = np.hstack((b1,-Q+np.dot(gQ,p)))

            # Constraint cut
            a2 = np.hstack((gS,0.,-1.,0.))
            A2 = np.vstack((A2,a2))
            b2 = np.hstack((b2,-S+np.dot(gS,self.x)))

            # Output
            if not quiet:
                print('{0:^8d}'.format(k), end=' ')
                print('{0:^10.2f}'.format(time.time()-t0), end=' ')
                print('{0:^12.5e}'.format(0.5*np.dot(p,H0*p)+np.dot(g0,p)+z1), end=' ')
                print('{0:^12.5e}'.format((1.-gamma)*t+z2), end=' ')
                print('{0:^12d}'.format(len(self.results)))
            
            # Update
            num_y += 1
            k += 1
            
            

        
        
