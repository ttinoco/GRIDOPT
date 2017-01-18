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
        import cvxpy as cpy
        from multiprocess import Pool

        # Parameters
        params = self.parameters
        num_sce = params['scenarios']
        num_procs = params['num_procs']
        maxiters = params['maxiters']
        maxtime = params['maxtime']
        period = params['period']
        z_inf = params['z_inf']
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

        # Header
        if not quiet:
            print('\nSAA Cutting-Plane Risk')
            print('-----------------------')
            print('{0:^8s}'.format('iter'), end=' ')
            print('{0:^10s}'.format('time(s)'), end=' ')
            print('{0:^12s}'.format('FL'), end= ' ')
            print('{0:^12s}'.format('GL'), end= ' ')
            print('{0:^12s}'.format('saved'))

        # Base problem
        p = cpy.Variable(num_p)
        t = cpy.Variable(1)
        z1 = cpy.Variable(1)
        z2 = cpy.Variable(1)
        obj = 0.5*cpy.quad_form(p,cpy.Constant(H0)) + p.T*cpy.Constant(g0) + z1
        constr = [p_min <= p, 
                  p <= p_max, 
                  t_min <= t, 
                  t <= t_max, 
                  -z_inf <= z1, 
                  -z_inf <= z2,
                  (1.-gamma)*t + z2 <= 0]

        # Init
        k = 0
        t1 = 0
        t0 = time.time()
        self.results = []

        # Loop
        while True:
            
            # Problem
            opt_prog = cpy.Problem(cpy.Minimize(obj),constr)
            
            # Solve
            opt_prog.solve(solver=cpy.ECOS,verbose=False,max_iters=1000)
                       
            # Results
            self.x = np.hstack((np.array(p.value).flatten(),t.value))
            assert(self.x.shape == (num_p+1,))
 
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
            chunk = int(np.ceil(float(num_sce)/float(num_procs)))
            Q_list,gQ_list = zip(*pool.map(lambda w: problem.eval_Q(self.x[:num_p],w),scenarios,chunksize=chunk))
            Q = sum(Q_list)/float(num_sce)
            gQ = sum(gQ_list)/float(num_sce)

            # Constr saa
            S_list = []
            gS_list = []
            for i in range(num_sce):
                qq = Q_list[i]-problemRA.Qmax-self.x[-1]
                S_list.append(np.maximum(qq,0.))
                if qq >= 0.:
                    gS_list.append(np.hstack((gQ_list[i],-1.)))
                else:
                    gS_list.append(np.hstack((op,0.)))
            S = sum(S_list)/float(num_sce)
            gS = sum(gS_list)/float(num_sce)

            # Obj cut
            constr.append(Q + (p-p.value).T*cpy.Constant(gQ) <= z1)

            # Constraint cut
            constr.append(S + cpy.vstack(p-p.value,t-t.value).T*cpy.Constant(gS) <= z2)

            # Output
            if not quiet:
                print('{0:^8d}'.format(k), end=' ')
                print('{0:^10.2f}'.format(time.time()-t0), end=' ')
                print('{0:^12.5e}'.format(0.5*np.dot(self.x[:num_p],H0*self.x[:num_p])+
                                          np.dot(g0,self.x[:num_p])+z1.value), end=' ')
                print('{0:^12.5e}'.format((1.-gamma)*self.x[-1]+z2.value), end=' ')
                print('{0:^12d}'.format(len(self.results)))
            
            # Update
            k += 1
            
            

        
        
