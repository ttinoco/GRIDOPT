#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from utils import ApplyFunc
from problem import TS_DCOPF
from scipy.sparse import csr_matrix,eye
from multiprocessing import Pool,cpu_count
from optalg.stoch_solver import StochGen_Problem

class TS_DCOPF_RiskAverse(StochGen_Problem):
    """"
    This class represents a problem of the form
    
    minimize(p,t)   varphi_0(p) + E_r[Q(p,r)]
    subject to      E_r[ (Q(p,r)-Qmax-t)_+ + (1-gamma)t ] <= 0 { or CVaR(Q(p,r)-Qmax,gamma) }
                    p_min <= p <= p_max
    
    where Q(p,r) is the optimal value of

    minimize(q,w,s,z)   varphi_1(q)
    subjcet to          G(p+q) + Rs - Aw = b
                        Jw - z = 0
                        p_min <= p+q <= p_max
                        z_min <= z <= z_max
                        0 <= s <= r.
    """

    # Parameters
    parameters = {'lam_max' : 1e3}   # max Lagrange multiplier
    
    def __init__(self,net,Qfac,gamma,samples):
        """
        Class constructor.
        
        Parameters
        ----------
        net : Network
        Qfac : float (> 1)
        gamma : float
        """

        # Parameters
        self.parameters = TS_DCOPF_RiskAverse.parameters.copy()
        
        # Save args
        self.Qfac = Qfac
        self.gamma = gamma
        
        # Regular problem
        self.ts_dcopf = TS_DCOPF(net)

        # Qmax
        p_ce = self.ts_dcopf.solve_approx(quiet=True)
        self.Qmax = Qfac*self.ts_dcopf.eval_EQ(p_ce,samples=samples)[0]

    def eval_FG(self,x,w,problem=None,debug=False):
        """
        Evaluates F, G and their subgradients at x
        for the given w.

        Parameters
        ----------
        x : (p,t)
        w : renewable powers

        Returns
        -------
        F : float
        gF : subgradient vector
        G : vector
        JG : csr_matrix of subgradients 
        """

        p = x[:-1]
        t = x[-1]

        gamma = self.gamma

        H0 = self.ts_dcopf.H0
        g0 = self.ts_dcopf.g0
        
        phi = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
        gphi = H0*p + g0
        Q,gQ = self.ts_dcopf.eval_Q(p,w,problem=problem)

        F =  phi+Q
        gF = np.hstack((gphi+gQ,0.))

        sigma = Q-self.Qmax-t

        G = np.array([np.maximum(sigma,0.) + (1.-gamma)*t])
        if sigma >= 0:
            JG = csr_matrix(np.hstack((gQ,-1. + (1.-gamma))),shape=(1,x.size))
        else:
            JG = csr_matrix(np.hstack((np.zeros(p.size),1.-gamma)),shape=(1,x.size))

        # Debug
        if debug:
            print Q,self.Qmax,t,sigma,self.ts_dcopf.get_prop_x(p)

        return F,gF,G,JG

    def eval_FG_approx(self,x):
        """
        Evaluates certainty-equivalent approximations
        of F and G and their derivaties.

        Parameters
        ----------
        x : (p,t)

        Returns
        -------
        F : float
        gF : gradient vector
        G : vector
        JG : Jacobian matrix
        """
        
        p = x[:-1]
        t = x[-1]
        Er = self.ts_dcopf.Er
        
        gamma = self.gamma

        H0 = self.ts_dcopf.H0
        g0 = self.ts_dcopf.g0
        
        phi = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
        gphi = H0*p + g0
        Q,gQ = self.ts_dcopf.eval_Q(p,Er)

        F =  phi+Q
        gF = np.hstack((gphi+gQ,0.))

        sigma = Q-self.Qmax-t
        if sigma > 100:
            C = 1.
            log_term = sigma
        else:
            C = np.exp(sigma)/(1.+np.exp(sigma))
            log_term = np.log(1. + np.exp(sigma))

        G = np.array([log_term + (1.-gamma)*t])
        JG = csr_matrix(np.hstack((C*gQ,-C + (1.-gamma))),shape=(1,x.size))

        return F,gF,G,JG

    def eval_EFG_sequential(self,x,samples=500):

        # Local vars
        p = x[:-1]
        t = x[-1]
        num_p = self.ts_dcopf.num_p
        num_w = self.ts_dcopf.num_w
        num_r = self.ts_dcopf.num_r

        # Seed
        np.random.seed()

        # Init
        F = 0.
        gF = np.zeros(x.size)
        G = np.zeros(1)
        JG = csr_matrix((1,x.size))
        
        # Second stage problem
        problem = self.ts_dcopf.get_problem_for_Q(p,np.ones(num_r))
        
        # Sampling loop
        for i in range(samples):
            
            r = self.sample_w()

            problem.u[num_p+num_w:num_p+num_w+num_r] = r # Important (update bound)
            
            F1,gF1,G1,JG1 = self.eval_FG(x,r,problem=problem)

            # Update
            F += (F1-F)/(i+1.)
            gF += (gF1-gF)/(i+1.)
            G += (G1-G)/(i+1.)
            JG = JG + (JG1-JG)/(i+1.)
                            
        return F,gF,G,JG
        
    def eval_EFG(self,x,samples=500,num_procs=None):

        if not num_procs:
            num_procs = cpu_count()
        pool = Pool(num_procs)
        num = int(np.ceil(float(samples)/float(num_procs)))
        results = zip(*pool.map(ApplyFunc,num_procs*[(self,'eval_EFG_sequential',x,num)]))
        return map(lambda vals: sum(map(lambda val: num*val/float(num*num_procs),vals)),results)
        
    def get_size_x(self):

        return self.ts_dcopf.num_p + 1

    def get_size_lam(self):

        return 1

    def get_prop_x(self,x):
        
        p = x[:-1]
        t = x[-1]
        
        return t#self.ts_dcopf.get_prop_x(p)
        
    def project_x(self,x):
        
        p = x[:-1]
        t = x[-1]

        return np.hstack((self.ts_dcopf.project_x(p),t))

    def project_lam(self,lam):

        lmax = self.parameters['lam_max']
        return np.maximum(np.minimum(lam,lmax),0.)

    def sample_w(self):

        return self.ts_dcopf.sample_w()
        
    def show(self):

        self.ts_dcopf.show()
        print 'Qmax : %.5e' %self.Qmax

    def solve_Lrelaxed_approx(self,lam,g_corr=None,J_corr=None,tol=1e-4,quiet=False):
        """
        Solves
        
        minimize(x)   F_approx + lam^TG_approx(x) + g^Tx + lam^TJx (slope correction)
        subject to    x in X

        Returns
        -------
        x : vector
        """

        # Local variables
        prob = self.ts_dcopf
        
        num_p = prob.num_p
        num_w = prob.num_w
        num_r = prob.num_r
        num_bus = prob.num_bus
        num_br = prob.num_br
        num_t = 1

        Ont = coo_matrix((num_bus,num_t))
        op = np.zeros(num_p)
        oz = np.zeros(num_br)

        Ip = eye(num_p,format='coo')
        Iz = eye(num_br,format='coo')

        Onp = coo_matrix((num_bus,num_p))
        Onz = coo_matrix((num_bus,num_br))        
        Ow = coo_matrix((num_w,num_w))
        Os = coo_matrix((num_r,num_r))
        Oz = coo_matrix((num_br,num_br))

        # Corrections
        if g_corr is None:
            pass # do something
        if J_corr is None:
            pass # do something
        
        # Problem construction
        A = bmat([[prob.G,Ont,prob.G,-prob.A,prob.R,None,None],
                  [Ip,None,Ip,None,None,-Ip,None],
                  [None,None,None,prob.J,None,None,-Iz]],format='coo')
        b = np.hstack((prob.b,op,oz))
    

