#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from utils import ApplyFunc
from types import MethodType
from numpy.linalg import norm
from problem import TS_DCOPF_Problem
from multiprocessing import Pool,cpu_count
from optalg.stoch_solver import StochGen_Problem
from optalg.opt_solver import OptProblem,OptSolverLCCP
from scipy.sparse import csr_matrix,eye,bmat,coo_matrix,tril

class TS_DCOPF_RA_Problem(StochGen_Problem):
    """"
    This class represents a problem of the form
    
    minimize(p,t)   varphi_0(p) + E_r[Q(p,r)]
    subject to      E_r[ (Q(p,r)-Qmax-t)_+ + (1-gamma)t ] <= 0 { or CVaR(Q(p,r)-Qmax,gamma) }
                    p_min <= p <= p_max
    
    where Q(p,r) is the optimal value of

    minimize(q,w,s,z)   varphi_1(q)
    subjcet to          G(p+q) + Rs - Aw = b
                        p_min <= p+q <= p_max
                        z_min <= Jw <= z_max
                        0 <= s <= r.
    """

    # Parameters
    parameters = {'lam_max' : 1e2,   # max Lagrange multiplier
                  'smax_param': 1e2, # softmax parameter
                  't_reg': 1e-3,
                  't_min': -0.1,
                  't_max': 0.}
    
    def __init__(self,net,Qfac,gamma,samples):
        """
        Class constructor.
        
        Parameters
        ----------
        net : Network
        Qfac : float (> 1)
        gamma : float
        samples : int
        """

        # Parameters
        self.parameters = TS_DCOPF_RA_Problem.parameters.copy()
        
        # Save args
        self.Qfac = Qfac    # factor for setting Qmax
        self.gamma = gamma  # parameter for CVaR (e.g. 0.95)
        
        # Regular problem
        self.ts_dcopf = TS_DCOPF_Problem(net)

        # Qnorm and Qmax
        p_ce,results = self.ts_dcopf.solve_approx(quiet=True)
        self.Qnorm = self.ts_dcopf.eval_EQ_parallel(p_ce,samples=samples)[0]
        self.Qmax = Qfac*self.Qnorm

        # Constants
        self.num_p = self.ts_dcopf.num_p
        self.num_w = self.ts_dcopf.num_w
        self.num_r = self.ts_dcopf.num_r
        self.num_bus = self.ts_dcopf.num_bus
        self.num_br = self.ts_dcopf.num_br
        self.temp_x = np.zeros(self.num_p+1) 
        self.JG_const = csr_matrix(np.hstack((np.zeros(self.num_p),1.-gamma)),shape=(1,self.temp_x.size))
        self.op = np.zeros(self.num_p)
        self.ow = np.zeros(self.num_w)
        self.os = np.zeros(self.num_r)
        self.oz = np.zeros(self.num_br)
        self.Ip = eye(self.num_p,format='coo')
        self.Iz = eye(self.num_br,format='coo')
        self.Ont = coo_matrix((self.num_bus,1))
        self.Ot = coo_matrix((1,1))
        self.Ow = coo_matrix((self.num_w,self.num_w))
        self.Os = coo_matrix((self.num_r,self.num_r))
        self.Op = coo_matrix((self.num_p,self.num_p))
        self.Oz = coo_matrix((self.num_br,self.num_br))
        self.ones_r = np.ones(self.num_r)
        self.ones_w = np.ones(self.num_w)
        
    def eval_VaR(self,p,t=0,iters=1000,tol=1e-4):
        """
        Evaluates (approximately) CVaR(Q(p,r)-Qmax,gamma).

        Parameters
        ----------
        p : vector
        t : float (initial estimate)

        Returns
        -------
        var : float
        """

        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r

        problem = self.ts_dcopf.get_problem_for_Q(p,self.ones_r)
        
        print 'var'
        print 'k     t'
        for k in range(iters):
          
            print '%5d    %.5e' %(k,t)   
            
            r = self.sample_w()
            
            problem.u[num_p+num_w:num_p+num_w+num_r] = r # Important (update bound)
            
            Q,gQ = self.ts_dcopf.eval_Q(p,r,problem=problem,tol=tol)
           
            if Q-self.Qmax-t >= 0.:
                g = 1. - 1./(1.-self.gamma)
            else:
                g = 1.

            t -= g/(k+1.)
        return t

    def eval_FG(self,x,w,problem=None,tol=1e-4,info=False):
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
        t_reg = self.parameters['t_reg']
        num_p = self.num_p
        num_x = num_p+1
        temp_x = self.temp_x

        H0 = self.ts_dcopf.H0
        g0 = self.ts_dcopf.g0
        
        phi0 = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
        gphi0 = H0*p + g0
        Q,gQ = self.ts_dcopf.eval_Q(p,w,problem=problem,tol=tol)

        F =  phi0+Q+0.5*t_reg*(t**2.)
        gF = np.hstack((gphi0+gQ,t_reg*t))
        
        ind = 1. if Q <= self.Qmax else 0.

        sigma = Q-self.Qmax-t
        
        G = np.array([np.maximum(sigma,0.) + (1.-gamma)*t])
        if sigma >= 0:
            temp_x[:-1] = gQ
            temp_x[-1] = -1. + (1.-gamma)
            JG = csr_matrix(temp_x,shape=(1,num_x))
        else:
            JG = self.JG_const

        if not info:
            return F,gF,G,JG
        else:
            return F,gF,G,JG,ind

    def eval_FG_approx(self,x,tol=1e-4):
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
        smax_param = self.parameters['smax_param']
        t_reg = self.parameters['t_reg']
        gamma = self.gamma
        
        H0 = self.ts_dcopf.H0
        g0 = self.ts_dcopf.g0
        
        phi0 = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
        gphi0 = H0*p + g0
        Q,gQ = self.ts_dcopf.eval_Q(p,Er,tol=tol)

        F =  phi0+Q+0.5*t_reg*(t**2.)
        gF = np.hstack((gphi0+gQ,t_reg*t))

        sigma = smax_param*(Q-self.Qmax-t)/self.Qnorm
        a = np.maximum(sigma,0.)
        C = np.exp(sigma-a)/(np.exp(-a)+np.exp(sigma-a))
        log_term = a + np.log(np.exp(-a) + np.exp(sigma-a))

        G = np.array([self.Qnorm*log_term/smax_param + (1.-gamma)*t])
        JG = csr_matrix(np.hstack((C*gQ,-C + 1.-gamma)),shape=(1,x.size))

        return F,gF,G,JG

    def eval_EFG_sequential(self,x,samples=500,seed=None,tol=1e-4,info=False):

        # Local vars
        p = x[:-1]
        t = x[-1]
        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r

        # Seed
        if seed is None:
            np.random.seed()
        else:
            np.random.seed(seed)

        # Init
        ind = 0.
        F = 0.
        gF = np.zeros(x.size)
        G = np.zeros(1)
        JG = csr_matrix((1,x.size))
        
        # Second stage problem
        problem = self.ts_dcopf.get_problem_for_Q(p,self.ones_r)
        
        # Sampling loop
        for i in range(samples):
            
            r = self.sample_w()
            
            problem.u[num_p+num_w:num_p+num_w+num_r] = r # Important (update bound)
            
            F1,gF1,G1,JG1,ind1 = self.eval_FG(x,r,problem=problem,tol=tol,info=True)

            # Update
            ind += (ind1-ind)/(i+1.)
            F += (F1-F)/(i+1.)
            gF += (gF1-gF)/(i+1.)
            G += (G1-G)/(i+1.)
            JG = JG + (JG1-JG)/(i+1.)
                 
        if not info:
            return F,gF,G,JG
        else:
            return F,gF,G,JG,ind
        
    def eval_EFG(self,x,samples=500,num_procs=None,tol=1e-4,info=False):

        if not num_procs:
            num_procs = cpu_count()
        pool = Pool(num_procs)
        num = int(np.ceil(float(samples)/float(num_procs)))
        results = zip(*pool.map(ApplyFunc,[(self,'eval_EFG_sequential',x,num,i,tol,info) for i in range(num_procs)],chunksize=1))
        pool.terminate()
        pool.join()
        return map(lambda vals: sum(map(lambda val: val/float(num_procs),vals)),results)
        
    def get_size_x(self):

        return self.num_p + 1

    def get_size_lam(self):

        return 1

    def get_prop_x(self,x):
        
        p = x[:-1]
        t = x[-1]
        
        return t#self.ts_dcopf.get_prop_x(p)
        
    def project_x(self,x):
        
        p = x[:-1]
        t = x[-1]
        t_max = self.parameters['t_max']*self.Qnorm
        t_min = self.parameters['t_min']*self.Qnorm

        return np.hstack((self.ts_dcopf.project_x(p),
                          np.maximum(np.minimum(t,t_max),t_min)))

    def project_lam(self,lam):

        lmax = self.parameters['lam_max']
        return np.maximum(np.minimum(lam,lmax),0.)

    def sample_w(self):

        return self.ts_dcopf.sample_w()

    def save_x_info(self,x,filename):
       
        self.ts_dcopf.save_x_info(x[:-1],filename)
 
    def show(self):

        self.ts_dcopf.show()
        print 'Qnorm      : %.5e' %self.Qnorm
        print 'Qmax       : %.5e' %self.Qmax
        print 'Qfac       : %.2f' %self.Qfac
        print 'gamma      : %.2f' %self.gamma
        print 'smax param : %.2e' %self.parameters['smax_param']
        print 'lmax       : %.2e' %self.parameters['lam_max']
        print 't_reg      : %.2e' %self.parameters['t_reg']
        print 't_min      : %.2e' %self.parameters['t_min']
        print 't_max      : %.2e' %self.parameters['t_max']

    def solve_Lrelaxed_approx(self,lam,g_corr=None,J_corr=None,tol=1e-4,quiet=False,init_data=None):
        """
        Solves
        
        minimize(x)   F_approx + lam^TG_approx(x) + g^Tx + lam^TJx (slope correction)
        subject to    x in X

        Returns
        -------
        x : vector
        """
        
        # Construct problem
        problem = self.construct_Lrelaxed_approx_problem(lam,g_corr=g_corr,J_corr=J_corr)

        # Warm start
        if init_data is not None:
            problem.x = init_data['x']
            problem.lam = init_data['lam']
            problem.mu = init_data['mu']
            problem.pi = init_data['pi']

        # Solve problem
        solver = OptSolverLCCP()
        solver.set_parameters({'quiet': quiet,
                               'tol': tol})
        try:
            solver.solve(problem)
        except Exception:
            raise
        finally:
            pass
        
        # Get results
        results = solver.get_results()
        x = results['x']
        lam = results['lam']
        mu = results['mu']
        pi = results['pi']
        
        # Check
        A = problem.A
        b = problem.b
        l = problem.l
        u = problem.u
        problem.eval(x)
        gphi = problem.gphi
        assert(norm(gphi-A.T*lam+mu-pi) < (1e-4)*(norm(gphi)+norm(lam)+norm(mu)+norm(pi)))
        assert(norm(mu*(u-x)) < (1e-4)*(norm(mu)+norm(u-x)))
        assert(norm(pi*(x-l)) < (1e-4)*(norm(pi)+norm(x-l)))
        assert(np.all(x < u + 1e-4))
        assert(np.all(x > l - 1e-4))
        assert(norm(A*x-b) < (1e-4)*norm(b))

        # Return
        return x[:self.ts_dcopf.num_p+1],results

    def construct_Lrelaxed_approx_problem(self,lam,g_corr=None,J_corr=None):

        # Local variables
        Qmax = self.Qmax
        Qnorm = self.Qnorm
        prob = self.ts_dcopf
        inf = prob.parameters['infinity']
        smax_param = self.parameters['smax_param']
        t_reg = self.parameters['t_reg']
        t_max = self.parameters['t_max']*Qnorm
        t_min = self.parameters['t_min']*Qnorm
        gamma = self.gamma
        lam = float(lam)
        
        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r
        num_bus = self.num_bus
        num_br = self.num_br

        H0 = prob.H0
        g0 = prob.g0
        H1 = prob.H1
        g1 = prob.g1
        
        op = self.op
        ow = self.ow
        os = self.os
        oz = self.oz
        
        Ip = self.Ip
        Iz = self.Iz
        
        Ont = self.Ont
        Ot = self.Ot
        Ow = self.Ow
        Os = self.Os
        Op = self.Op
        Oz = self.Oz

        # Corrections
        if g_corr is None:
            g_corr = np.zeros(num_p+1)
        if J_corr is None:
            J_corr = np.zeros(num_p+1)
        else:
            J_corr = J_corr.toarray()[0,:]
        eta_p = g_corr[:-1]
        eta_t = g_corr[-1]
        nu_p = J_corr[:-1]
        nu_t = J_corr[-1]        
        
        # Problem construction
        A = bmat([[prob.G,Ont,prob.G,-prob.A,prob.R,None,None],
                  [Ip,None,Ip,None,None,-Ip,None],
                  [None,None,None,prob.J,None,None,-Iz]],format='coo')
        b = np.hstack((prob.b,op,oz))
        l = np.hstack((prob.p_min,           # p
                       t_min,                # t
                       -prob.p_max+prob.p_min, # q
                       -inf*self.ones_w,     # theta
                       self.os,              # s
                       prob.p_min,           # y
                       prob.z_min))          # z
        u = np.hstack((prob.p_max,           # p
                       t_max,                # t
                       prob.p_max-prob.p_min,  # q
                       inf*self.ones_w,      # theta
                       prob.Er,              # s
                       prob.p_max,           # y
                       prob.z_max))          # z
                
        def eval(cls,x):
            
            # Extract components
            offset = 0
            p = x[offset:offset+num_p]
            offset += num_p
            
            t = x[offset]
            offset += 1
            
            q = x[offset:offset+num_p]
            offset += num_p
            
            w = x[offset:offset+num_w]
            offset += num_w
            
            s = x[offset:offset+num_r]
            offset += num_r

            y = x[offset:offset+num_p]
            offset += num_p

            z = x[offset:]
            assert(z.size == num_br)

            # Eval partial functions
            phi0 = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
            gphi0 = H0*p + g0

            phi1 = 0.5*np.dot(q,H1*q)+np.dot(g1,q)
            gphi1 = H1*q + g1
            
            beta = smax_param*(phi1-Qmax-t)/Qnorm
            a = np.maximum(beta,0.)
            ebma = np.exp(beta-a)
            ebm2a = np.exp(beta-2*a)
            ema = np.exp(-a)
            C1 = ebma/(ema+ebma)
            C2 = smax_param*ebm2a/(Qnorm*(ema*ema+2*ebm2a+ebma*ebma))
            log_term = a + np.log(ema+ebma)
            
            # Value
            cls.phi = (phi0 + 
                       phi1 +
                       0.5*t_reg*(t**2.)+
                       lam*Qnorm*log_term/smax_param + lam*(1-gamma)*t + 
                       np.dot(eta_p+lam*nu_p,p) + 
                       (eta_t+lam*nu_t)*t)
            
            # Gradient
            cls.gphi = np.hstack((gphi0 + eta_p + lam*nu_p, # p
                                  t_reg*t + lam*(-C1 + (1.-gamma)) + eta_t + lam*nu_t, # t
                                  (1.+lam*C1)*gphi1, # q
                                  ow,                     # theta
                                  os,                     # s
                                  op,                     # y
                                  oz))                    # z
            
            # Hessian (lower triangular)
            H = (1.+lam*C1)*H1 + tril(lam*C2*np.outer(gphi1,gphi1))
            g = gphi1.reshape((q.size,1))
            cls.Hphi = bmat([[H0,None,None,None,None,None,None],             # p
                             [None,t_reg+lam*C2,-lam*C2*g.T,None,None,None,None],  # t
                             [None,-lam*C2*g,H,None,None,None,None],         # q
                             [None,None,None,Ow,None,None,None],         # theta
                             [None,None,None,None,Os,None,None],         # s
                             [None,None,None,None,None,Op,None],         # y
                             [None,None,None,None,None,None,Oz]],        # z
                            format='coo')
            
        problem = OptProblem()
        problem.A = A
        problem.b = b
        problem.u = u
        problem.l = l
        problem.eval = MethodType(eval,problem)
        
        return problem


