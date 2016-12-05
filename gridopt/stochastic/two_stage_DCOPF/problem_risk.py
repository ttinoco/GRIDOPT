#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from scipy.sparse import csr_matrix
from .problem import TS_DCOPF_Problem
from optalg.stoch_solver import StochProblemC

class TS_DCOPF_RA_Problem(StochProblemC):
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
    parameters = {'lam_max' : 1e2,      # max Lagrange multiplier
                  't_reg': 1e-8,
                  't_min': -0.1,
                  't_max': 0.,
                  'Qfac': 0.8,          # factor for setting Qmax
                  'gamma': 0.95,        # parameter for CVaR (e.g. 0.95)
                  'num_samples' : 1000,
                  'num_procs': 1}
    
    def __init__(self,net,parameters={}):
        """
        Class constructor.
        
        Parameters
        ----------
        net : Network
        parameters : dict
        """

        # Parameters
        self.parameters = TS_DCOPF_RA_Problem.parameters.copy()
        self.set_parameters(parameters)
        
        # Local vars
        Qfac = self.parameters['Qfac']
        gamma = self.parameters['gamma']
 
        # Regular problem
        self.ts_dcopf = TS_DCOPF_Problem(net,self.parameters)

        # Qref and Qmax
        p_ce,gF_ce = self.ts_dcopf.solve_approx()
        self.Qref = self.ts_dcopf.eval_EQ(p_ce)[0]
        self.Fref = 0.5*np.dot(p_ce,self.ts_dcopf.H0*p_ce)+np.dot(self.ts_dcopf.g0,p_ce)+self.Qref
        self.Qmax = Qfac*self.Qref

        # Constants
        self.num_p = self.ts_dcopf.num_p
        self.num_w = self.ts_dcopf.num_w
        self.num_r = self.ts_dcopf.num_r
        self.temp_x = np.zeros(self.num_p+1) 
        self.JG_const = csr_matrix(np.hstack((np.zeros(self.num_p),1.-gamma)),shape=(1,self.temp_x.size))
        
        # Opt prog
        self.clear()

    def clear(self):

        self.ts_dcopf.clear()
        self.opt_prog = None
        self.opt_prog_data = None

    def set_parameters(self,params):
        """
        Sets problem parameters.
        
        Parameters
        ----------
        params : dict
        """
        
        for key,value in list(params.items()):
            if key in self.parameters:
                self.parameters[key] = value
        
    def eval_FG(self,x,w,info=False):
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

        gamma = self.parameters['gamma']
        temp_x = self.temp_x

        H0 = self.ts_dcopf.H0
        g0 = self.ts_dcopf.g0
        
        phi0 = 0.5*np.dot(p,H0*p)+np.dot(g0,p)
        gphi0 = H0*p + g0
        Q,gQ = self.ts_dcopf.eval_Q(p,w)

        F =  phi0+Q
        gF = np.hstack((gphi0+gQ,0.))
        
        ind = 1. if Q <= self.Qmax else 0.

        sigma = Q-self.Qmax-t
        
        G = np.array([np.maximum(sigma,0.) + (1.-gamma)*t])
        if sigma >= 0:
            temp_x[:-1] = gQ
            temp_x[-1] = -1. + (1.-gamma)
            JG = csr_matrix(temp_x,shape=(1,temp_x.size))
        else:
            JG = self.JG_const

        if not info:
            return F,gF,G,JG
        else:
            return F,gF,G,JG,ind
        
    def eval_EFG_sequential(self,x,num_samples=500,seed=None,info=False):
        
        # Local vars
        p = x[:-1]
        t = x[-1]
 
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
                
        # Sampling loop
        for i in range(num_samples):
            
            r = self.sample_w()
                        
            F1,gF1,G1,JG1,ind1 = self.eval_FG(x,r,info=True)

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
        
    def eval_EFG(self,x,info=False):

        from multiprocess import Pool

        num_procs = self.parameters['num_procs']
        num_samples = self.parameters['num_samples']
        pool = Pool(num_procs)
        num = int(np.ceil(float(num_samples)/float(num_procs)))
        results = list(zip(*map(lambda i: self.eval_EFG_sequential(x,num,i,info),range(num_procs))))
        pool.terminate()
        pool.join()
        if not info:
            assert(len(results) == 4)
        else:
            assert(len(results) == 5)
        assert(all([len(vals) == num_procs for vals in results]))
        return [sum(vals)/float(num_procs) for vals in results]
        
    def get_size_x(self):

        return self.num_p + 1

    def get_init_x(self):

        x0,gF_approx,JG_approx,results = self.solve_Lrelaxed_approx(np.zeros(self.get_size_lam()),quiet=True)
        x0[-1] = self.parameters['t_min']*self.Qref
        return x0

    def get_size_lam(self):

        return 1

    def get_prop_x(self,x):
        
        p = x[:-1]
        t = x[-1]
        
        return t #self.ts_dcopf.get_prop_x(p)
        
    def project_x(self,x):
        
        p = x[:-1]
        t = x[-1]
        t_max = self.parameters['t_max']*self.Qref
        t_min = self.parameters['t_min']*self.Qref

        return np.hstack((self.ts_dcopf.project_x(p),
                          np.maximum(np.minimum(t,t_max),t_min)))

    def project_lam(self,lam):

        lmax = self.parameters['lam_max']
        return np.maximum(np.minimum(lam,lmax),0.)

    def sample_w(self):

        return self.ts_dcopf.sample_w()
 
    def show(self):

        self.ts_dcopf.show()

        print('Fref        : %.5e' %self.Fref)
        print('Qref        : %.5e' %self.Qref)
        print('Qmax        : %.5e' %self.Qmax)
        print('Qfac        : %.2f' %self.parameters['Qfac'])
        print('gamma       : %.2f' %self.parameters['gamma'])
        print('lmax        : %.2e' %self.parameters['lam_max'])
        print('t_reg       : %.2e' %self.parameters['t_reg'])
        print('t_min       : %.2e' %self.parameters['t_min'])
        print('t_max       : %.2e' %self.parameters['t_max'])
        print('num_samples : %d' %self.parameters['num_samples'])
        print('num procs   : %d' %self.parameters['num_procs'])

    def solve_Lrelaxed_approx(self,lam,g_corr=None,J_corr=None,quiet=True,init_data=None):
        """
        Solves
        
        minimize(x)   F_approx + lam^TG_approx(x) + g^Tx + lam^TJx (slope correction)
        subject to    x in X

        Parameters
        ----------

        Returns
        -------
        x : vector
        """

        # Local vars
        t_reg = self.parameters['t_reg']
        gamma = self.parameters['gamma']
        prob = self.ts_dcopf

        # Multiplier
        lam = float(lam)

        # Corrections
        if g_corr is None:
            g_corr = np.zeros(self.num_p+1)
        if J_corr is None:
            J_corr = np.zeros(self.num_p+1)
        else:
            J_corr = J_corr.toarray()[0,:]
        eta_p = g_corr[:-1]
        eta_t = g_corr[-1]
        nu_p = J_corr[:-1]
        nu_t = J_corr[-1]

        # Opt program
        opt_prog,opt_prog_data = self.get_Lrelaxed_approx_program()
        opt_prog_data['lam'].value = lam
        opt_prog_data['eta_p'].value = eta_p
        opt_prog_data['eta_t'].value = eta_t
        opt_prog_data['nu_p'].value = nu_p
        opt_prog_data['nu_t'].value = nu_t

        # Solve
        opt_prog.solve()
           
        # Results
        p = np.array(opt_prog_data['p'].value).flatten()
        t = float(np.array(opt_prog_data['t'].value).flatten())
        q = np.array(opt_prog_data['q'].value).flatten()
        
        # Gradients
        Q = 0.5*np.dot(q,prob.H1*q)+np.dot(prob.g1,q)
        gQ = -(prob.H1*q+prob.g1)
        gF_approx = np.hstack((prob.H0*p+prob.g0+gQ,t_reg*t))
        if Q-self.Qmax-t >= 0.:
            JG_approx = csr_matrix(np.hstack((gQ,-1. + 1.-gamma)),shape=(1,p.size+1))
        else:
            JG_approx = csr_matrix(np.hstack((np.zeros(p.size),1.-gamma)),shape=(1,p.size+1))
        return np.hstack((p,t)),gF_approx,JG_approx,None
        
    def get_Lrelaxed_approx_program(self):

        # Construction
        if self.opt_prog is None:

            # Import
            import cvxpy as cpy
            
            # Local vars
            Qmax = self.Qmax
            Qref = self.Qref
            t_reg = self.parameters['t_reg']
            t_max = self.parameters['t_max']*Qref
            t_min = self.parameters['t_min']*Qref
            gamma = self.parameters['gamma']
            prob = self.ts_dcopf
            num_p = self.num_p
            num_w = self.num_w
            num_r = self.num_r
            
            H0 = prob.H0
            g0 = prob.g0
            H1 = prob.H1
            g1 = prob.g1
            
            G = prob.G
            A = prob.A
            b = prob.b
            R = prob.R
            J = prob.J
            
            p_min = prob.p_min
            p_max = prob.p_max
            z_min = prob.z_min
            z_max = prob.z_max

            Er = prob.Er

            # Variables
            p = cpy.Variable(num_p)
            t = cpy.Variable(1)
            q = cpy.Variable(num_p)
            w = cpy.Variable(num_w)
            s = cpy.Variable(num_r)
                
            # Parameters
            lam_par = cpy.Parameter(sign="positive")
            eta_p_par = cpy.Parameter(num_p)
            eta_t_par = cpy.Parameter()
            nu_p_par = cpy.Parameter(num_p)
            nu_t_par = cpy.Parameter()

            obj = (0.5*cpy.quad_form(p,cpy.Constant(H0)) + p.T*cpy.Constant(g0) +
                   0.5*cpy.quad_form(q,cpy.Constant(H1)) + q.T*cpy.Constant(g1) +
                   0.5*t_reg*cpy.square(t) +
                   p.T*eta_p_par +
                   t*eta_t_par +
                   lam_par*((1.-gamma)*t + 
                            cpy.max_elemwise(0.5*cpy.quad_form(q,cpy.Constant(H1)) + q.T*cpy.Constant(g1)-Qmax-t,0) +
                            p.T*nu_p_par +
                            t*nu_t_par))
            constr = [p_min <= p, 
                      p <= p_max,
                      t_min <= t,
                      t <= t_max,
                      p_min <= p + q,
                      p + q <= p_max,
                      0 <= s,
                      s <= prob.Er,
                      z_min <= cpy.Constant(J)*w,
                      cpy.Constant(J)*w <= z_max,
                      cpy.Constant(G)*(p+q) + cpy.Constant(R)*s - cpy.Constant(A)*w == b]
            self.opt_prog = cpy.Problem(cpy.Minimize(obj),constr)
            self.opt_prog_data = {'lam': lam_par,
                                  'eta_p':eta_p_par,
                                  'eta_t':eta_t_par,
                                  'nu_p':nu_p_par,
                                  'nu_t':nu_t_par,
                                  'p':p,
                                  't':t,
                                  'q':q}
        
        # Return
        return self.opt_prog,self.opt_prog_data
