#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet as pf
import numpy as np
from numpy.linalg import norm
from scipy.sparse import triu,coo_matrix
from optalg.stoch_solver import StochProblem
            
class TS_DCOPF_Problem(StochProblem):
    """"
    This class represents a problem of the form
    
    minimize(p)   varphi_0(p) + E_r[Q(p,r)]
    subject to    p_min <= p <= p_max
    
    where Q(p,r) is the optimal value of

    minimize(q,w,s,z)   varphi_1(q)
    subjcet to          G(p+q) + Rs - Aw = b
                        Jw - z = 0
                        p_min <= p+q <= p_max
                        z_min <= z <= z_max
                        0 <= s <= r.
    """

    # Parameters
    parameters = {'cost_factor' : 1e2,   # factor for determining gen adjustment cost
                  'infinity' : 1e4,      # infinity
                  'flow_factor' : 1.,    # factor for relaxing thermal limits
                  'num_samples' : 1000,  # number of samples
                  'num_procs': 1}

    def __init__(self,net,parameters={}):
        """
        Class constructor.
        
        Parameters
        ----------
        net : PFNET Network
        parameters : dict
        """

        # Parameters
        self.parameters = TS_DCOPF_Problem.parameters.copy()
        self.set_parameters(parameters)
        
        # Save info
        self.total_load = sum([l.P for l in net.loads])
        self.uncertainty = 100.*sum([g.P_std for g in net.var_generators])/sum([g.P_max for g in net.var_generators]) # % of capacity
        self.corr_value = net.vargen_corr_value    # correlation value  ([0,1])
        self.corr_radius = net.vargen_corr_radius  # correlation radius (# of branches)
                
        # Branch flow limits
        for br in net.branches:
            if br.ratingA == 0.:
                br.ratingA = self.parameters['infinity']
            else:
                br.ratingA *= self.parameters['flow_factor']
        
        # Counters
        num_w = net.num_buses-net.get_num_slack_buses() # voltage angles
        num_p = net.get_num_P_adjust_gens()             # adjustable generators
        num_r = net.num_var_generators                  # renewable generators

        # Variables
        net.clear_flags()
        net.set_flags(pf.OBJ_BUS,
                      pf.FLAG_VARS,
                      pf.BUS_PROP_NOT_SLACK,
                      pf.BUS_VAR_VANG)
        net.set_flags(pf.OBJ_GEN,
                      pf.FLAG_VARS,
                      pf.GEN_PROP_P_ADJUST,
                      pf.GEN_VAR_P)
        net.set_flags(pf.OBJ_VARGEN,
                      pf.FLAG_VARS,
                      pf.VARGEN_PROP_ANY,
                      pf.VARGEN_VAR_P)
        assert(net.num_vars == num_w+num_p+num_r)
        
        # Values
        x = net.get_var_values()

        # Projections
        Pw = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)
        Pp = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P)
        Pr = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)
        assert(Pw.shape == (num_w,net.num_vars))
        assert(Pp.shape == (num_p,net.num_vars))
        assert(Pr.shape == (num_r,net.num_vars))

        # Power flow equations
        pf_eq = pf.Constraint(pf.CONSTR_TYPE_DCPF,net)
        pf_eq.analyze()
        pf_eq.eval(x)
        A = pf_eq.A.copy()
        b = pf_eq.b.copy()

        # Branch flow limits
        fl_lim = pf.Constraint(pf.CONSTR_TYPE_DC_FLOW_LIM,net)
        fl_lim.analyze()
        fl_lim.eval(x)
        G = fl_lim.G.copy()
        hl = fl_lim.l.copy()
        hu = fl_lim.u.copy()
        assert(np.all(hl < hu))
        
        # Generation cost
        cost = pf.Function(pf.FUNC_TYPE_GEN_COST,1.,net)
        cost.analyze()
        cost.eval(x)
        H = (cost.Hphi + cost.Hphi.T - triu(cost.Hphi))/net.base_power # symmetric
        g = cost.gphi/net.base_power - H*x
        assert(norm((H-H.T).data) < 1e-12*norm(H.data))
        
        # Bounds
        l = net.get_var_values(pf.LOWER_LIMITS)
        u = net.get_var_values(pf.UPPER_LIMITS)
        assert(np.all(Pw*l < Pw*u))
        assert(np.all(Pp*l < Pp*u))
        assert(np.all(Pr*l < Pr*u))
        
        # Save problem data
        self.num_w = num_w
        self.num_p = num_p
        self.num_r = num_r
        self.p_max = Pp*u
        self.p_min = Pp*l
        self.r_max = Pr*u
        self.r_base = Pr*x
        self.z_max = hu
        self.z_min = hl 
        self.H0 = Pp*H*Pp.T
        self.g0 = Pp*g
        self.H1 = self.H0*self.parameters['cost_factor']
        self.g1 = np.zeros(num_p)
        self.G = A*Pp.T
        self.R = A*Pr.T
        self.A = -A*Pw.T
        self.J = G*Pw.T
        self.b = b
        self.Pp = Pp
        self.Pw = Pw
        self.Pr = Pr

        # Renewable covariance
        from scikits.sparse.cholmod import cholesky
        r_cov = Pr*net.create_vargen_P_sigma(self.corr_radius,self.corr_value)*Pr.T
        self.r_cov = (r_cov+r_cov.T-triu(r_cov)).tocsc()
        assert(self.r_cov.shape == (self.num_r,self.num_r))
        factor = cholesky(self.r_cov)
        L,D = factor.L_D()
        P = factor.P()
        PT = coo_matrix((np.ones(P.size),(P,np.arange(P.size))),shape=D.shape)
        P = P.T
        D = D.tocoo()
        Dh = coo_matrix((np.sqrt(D.data),(D.row,D.col)),shape=D.shape)
        self.L = PT*L*Dh
        assert(np.all(D.row == D.col))
        assert(np.all(Dh.row == Dh.col))
        assert(np.all(D.data > 0))
        assert(np.all(Dh.data > 0))
        for i in range(10):
            z = np.random.randn(self.num_r)
            assert(norm(self.r_cov*z-self.L*self.L.T*z) < 1e-10)

        # Average renewables
        self.Er = np.zeros(self.num_r)
        for i in range(self.parameters['num_samples']):
            self.Er += (self.sample_w()-self.Er)/(i+1)
        
        # Checks
        assert(np.all(self.p_min == 0))
        assert(np.all(0 < self.r_max))
        assert(np.all(self.z_min < self.z_max))
        assert(np.all(self.p_min < self.p_max))
        assert(np.all(cost.Hphi.row == cost.Hphi.col))
        assert(np.all(cost.Hphi.data > 0))
        assert(norm(self.A.T*np.ones(net.num_buses)) < (1e-10)*np.sqrt(float(net.num_buses)))

        # Opt prog
        self.clear()
        
    def clear(self):

        self.opt_progQ = None
        self.opt_progQ_data = None

        self.opt_progF = None
        self.opt_progF_data = None
        
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
        
    def eval_EQ_sequential(self,p,num_samples=500,seed=None):
        """
        Evaluates E[Q(p,r)] and its gradient. 

        Parameters
        ----------
        p : generator powers
        samples : number of samples
        seed : integer
        """
        
        # Local vars
        Q = 0.
        gQ = np.zeros(self.num_p)

        # Seed
        if seed is None:
            np.random.seed()
        else:
            np.random.seed(seed)
                
        # Sampling loop
        for i in range(num_samples):
            
            r = self.sample_w()
            
            q,gq = self.eval_Q(p,r)

            # Update
            Q += (q-Q)/(i+1.)
            gQ += (gq-gQ)/(i+1.)
                            
        return Q,gQ

    def eval_EQ(self,p):
        """
        Evaluates E[Q(p,r)] and its gradient in parallel. 

        Parameters
        ----------
        p : generator powers
        """
       
        from multiprocess import Pool
 
        self.clear()
        num_procs = self.parameters['num_procs']
        num_samples = self.parameters['num_samples']
        pool = Pool(num_procs)
        num = int(np.ceil(float(num_samples)/float(num_procs)))
        results = list(zip(*pool.map(lambda i: self.eval_EQ_sequential(p,num,i),range(num_procs))))
        pool.terminate()
        pool.join()
        assert(len(results) == 2)
        assert(all([len(vals) == num_procs for vals in results]))
        return [sum(vals)/float(num_procs) for vals in results]
        
    def eval_Q(self,p,r):
        """
        Evaluates Q(p,r) and its gradient.

        Parameters
        ----------
        p : generator powers
        r : renewable powers

        Returns
        -------
        Q : value
        gQ : gradient
        """
        
        # Check
        assert(np.all(r > 0))
        
        # Program
        opt_progQ,opt_progQ_data = self.get_Q_program()
        opt_progQ_data['p'].value = p
        opt_progQ_data['r'].value = r
                           
        # Solve
        try: 
            opt_progQ.solve(max_iters=1000)
        except Exception:
            import cvxpy as cpy
            opt_progQ.solve(solver=cpy.CVXOPT,max_iters=1000)

        # Results
        q = np.array(opt_progQ_data['q'].value).flatten()
        
        # Objective
        Q = 0.5*np.dot(q,self.H1*q)+np.dot(self.g1,q)
        
        # Gradient
        gQ = -(self.H1*q+self.g1) # See ECC paper
        
        # Return
        return Q,gQ

    def get_size_x(self):
        """
        Gets size of vector of first stage variables.
        
        Returns
        -------
        size : int
        """

        return self.num_p

    def get_prop_x(self,x):
        """
        Gets some property of x that is useful for monitor 
        progress.

        Parameters
        ----------
        x : generator powers

        Returns
        -------
        prop : float
        """
        
        return np.average(x/self.p_max)
        
    def get_Q_program(self):
        """
        Constructs second-stage quadratic problem
        in the form

        minimize(x)   (1/2)x^THx + g^Tx
        subject to    Ax = b
                      l <= x <= u,

        where x = (q,w,s,z).

        Returns
        -------
        opt_progQ
        opt_progQ_data
        """
                
        # Construction
        if self.opt_progQ is None:
        
            # Import
            import cvxpy as cpy
    
            # Opt prog
            q = cpy.Variable(self.num_p)
            w = cpy.Variable(self.num_w)
            s = cpy.Variable(self.num_r)
            ppar = cpy.Parameter(self.num_p)
            rpar = cpy.Parameter(self.num_r)
            
            obj = 0.5*cpy.quad_form(q,cpy.Constant(self.H1)) + q.T*cpy.Constant(self.g1)
            
            constr = [self.p_min <= ppar + q,
                      ppar + q <= self.p_max,
                      0 <= s,
                      s <= rpar,
                      self.z_min <= cpy.Constant(self.J)*w,
                      cpy.Constant(self.J)*w <= self.z_max,
                      cpy.Constant(self.G)*(ppar+q) + cpy.Constant(self.R)*s - cpy.Constant(self.A)*w == self.b]
            
            self.opt_progQ = cpy.Problem(cpy.Minimize(obj),constr)
            self.opt_progQ_data = {'q': q,
                                   'p': ppar,
                                   'r': rpar}
        
        # Return
        return self.opt_progQ,self.opt_progQ_data
        
    def eval_F(self,x,w):
        """
        Evaluates objective function for a given
        realization of uncertainty.

        Parameters
        ----------
        x : generator powers
        w : renewable powers

        Returns
        -------
        F : float
        gF : vector
        """
        
        phi = 0.5*np.dot(x,self.H0*x)+np.dot(self.g0,x)
        gphi = self.H0*x + self.g0
        Q,gQ = self.eval_Q(x,w)

        return phi+Q,gphi+gQ

    def eval_EF(self,x):
        """
        Evaluates expected objective function.

        Parameters
        ----------
        x : generator powers
        
        Returns
        -------
        val : float
        """
        
        phi = 0.5*np.dot(x,self.H0*x)+np.dot(self.g0,x)
        gphi = self.H0*x + self.g0
        Q,gQ = self.eval_EQ(x)

        return phi+Q,gphi+gQ

    def project_x(self,x):
        """
        Projects generator powers on the set
        defined by generator limits.
        
        Parameters
        ----------
        x : generator powers

        Returns
        -------
        y : vector
        """

        return np.maximum(np.minimum(x,self.p_max),self.p_min)

    def sample_w(self):
        """
        Samples renewable powers.

        Returns
        -------
        r : vector
        """

        return np.minimum(np.maximum(self.r_base+self.L*np.random.randn(self.num_r),1e-3),self.r_max)

    def show(self):
        """
        Shows problem information.
        """

        Ctot = np.sum(self.r_max)
        Btot = np.sum(self.r_base)
        
        print('Stochastic Two-Stage DCOPF')
        print('--------------------------')
        print('gens             : %d' %self.num_p)
        print('angles           : %d' %self.num_w)
        print('vargens          : %d' %self.num_r)
        print('penetration cap  : %.2f (%% of load)' %(100.*Ctot/self.total_load))
        print('penetration base : %.2f (%% of load)' %(100.*Btot/self.total_load))
        print('penetration std  : %.2f (%% of local cap)' %self.uncertainty)
        print('correlation rad  : %d (edges)' %self.corr_radius)
        print('correlation val  : %.2f (unitless)' %self.corr_value)
        print('num samples      : %d' %self.parameters['num_samples'])
        print('num procs        : %d' %self.parameters['num_procs'])

    def get_F_program(self):
        
        # Construction
        if self.opt_progF is None:

            # Import
            import cvxpy as cpy
            
            # Opt prog
            p = cpy.Variable(self.num_p)
            q = cpy.Variable(self.num_p)
            w = cpy.Variable(self.num_w)
            s = cpy.Variable(self.num_r)
            gpar = cpy.Parameter(self.num_p)
            
            obj = (0.5*cpy.quad_form(p,cpy.Constant(self.H0)) + p.T*cpy.Constant(self.g0) +
                   p.T*gpar + 
                   0.5*cpy.quad_form(q,cpy.Constant(self.H1)) + q.T*cpy.Constant(self.g1))                   
            
            constr = [self.p_min <= p,
                      p <= self.p_max,
                      self.p_min <= p + q,
                      p + q <= self.p_max,
                      0 <= s,
                      s <= self.Er,
                      self.z_min <= cpy.Constant(self.J)*w,
                      cpy.Constant(self.J)*w <= self.z_max,
                      cpy.Constant(self.G)*(p+q) + cpy.Constant(self.R)*s - cpy.Constant(self.A)*w == self.b]
            
            self.opt_progF = cpy.Problem(cpy.Minimize(obj),constr)
            self.opt_progF_data = {'p': p,
                                   'q': q,
                                   'g': gpar}
        
        # Return
        return self.opt_progF,self.opt_progF_data

    def solve_approx(self,g_corr=None):
        """
        Solves certainty equivalent problem
        with slope correction.

        Parameters
        ----------
        g_corr : slope correction
        
        Returns
        -------
        x : vector
        gF_approx : gradient
        """
                
        # Slope correction
        if g_corr is None:
            g_corr = np.zeros(self.num_p)
            
        # Program
        opt_progF,opt_progF_data = self.get_F_program()
        opt_progF_data['g'].value = g_corr
                           
        # Solve
        opt_progF.solve()
        
        # Results
        p = np.array(opt_progF_data['p'].value).flatten()
        q = np.array(opt_progF_data['q'].value).flatten()
                
        # Gradient
        gF_approx = (self.H0*p+self.g0)-(self.H1*q+self.g1)
        
        # Returns
        return p,gF_approx
