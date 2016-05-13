#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import csv
import pfnet as pf
import numpy as np
from utils import ApplyFunc
from numpy.linalg import norm
from multiprocessing import Pool,cpu_count
from optalg.lin_solver import new_linsolver
from scipy.sparse.linalg import LinearOperator
from optalg.opt_solver.opt_solver_error import *
from optalg.stoch_solver import StochObj_Problem
from optalg.opt_solver import OptSolverIQP,QuadProblem
from scipy.sparse import triu,bmat,coo_matrix,eye,spdiags
            
class TS_DCOPF_Problem(StochObj_Problem):
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
                  'infinity' : 1e3,      # infinity
                  'flow_factor' : 1.,    # factor for relaxing thermal limits
                  'num_samples' : 2000}  # number of samples

    def __init__(self,net):
        """
        Class constructor.
        
        Parameters
        ----------
        net : PFNET Network
        """

        # Parameters
        self.parameters = TS_DCOPF_Problem.parameters.copy()

        # Save info
        self.total_load = sum([l.P for l in net.loads])
        self.uncertainty = 100.*sum([g.P_std for g in net.var_generators])/sum([g.P_max for g in net.var_generators]) # % of capacity
        self.corr_value = net.vargen_corr_value    # correlation value  ([0,1])
        self.corr_radius = net.vargen_corr_radius  # correlation radius (# of branches)
                
        # Generator limits
        for gen in net.generators:
            gen.P_min = 0.
            gen.P_max = np.maximum(gen.P_max,0.)
            assert(gen.P_min <= gen.P_max)

        # Branch flow limits
        for br in net.branches:
            if br.ratingA == 0.:
                br.ratingA = self.parameters['infinity']
            else:
                br.ratingA *= self.parameters['flow_factor']
        
        # Counters
        num_w = net.num_buses-net.get_num_slack_buses() # voltage angles
        num_p = net.get_num_P_adjust_gens()             # adjustable generators
        num_r = net.num_vargens                         # renewable generators
        num_bus = net.num_buses                         # buses
        num_br = net.num_branches                       # branches

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
        self.num_bus = num_bus
        self.num_br = num_br
        self.p_max = Pp*u
        self.p_min = Pp*l
        self.w_max = self.parameters['infinity']*np.ones(self.num_w)
        self.w_min = -self.parameters['infinity']*np.ones(self.num_w)
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
        self.Er = 0
        for i in range(self.parameters['num_samples']):
            self.Er += (self.sample_w()-self.Er)/(i+1)
        
        # Checks
        assert(np.all(self.p_min == 0))
        assert(np.all(0 < self.r_max))
        assert(np.all(self.z_min < self.z_max))
        assert(np.all(self.p_min < self.p_max))
        assert(np.all(cost.Hphi.row == cost.Hphi.col))
        assert(np.all(cost.Hphi.data > 0))
        assert(norm(self.A.T*np.ones(self.num_bus)) < (1e-10)*np.sqrt(self.num_bus*1.))
        
    def eval_EQ(self,p,samples=500,seed=None,tol=1e-4,quiet=True):
        """
        Evaluates E[Q(p,r)] and its gradient. 

        Parameters
        ----------
        p : generator powers
        samples : number of samples
        seed : integer
        tol : evaluation tolerance
        quiet : flag
        """
        
        # Local vars
        Q = 0.
        gQ = np.zeros(self.num_p)

        # Seed
        if seed is None:
            np.random.seed()
        else:
            np.random.seed(seed)
        
        # Problem
        problem = self.get_problem_for_Q(p,np.ones(self.num_r))
        
        # Sampling loop
        for i in range(samples):
            
            r = self.sample_w()

            problem.u[self.num_p+self.num_w:self.num_p+self.num_w+self.num_r] = r # Important (update bound)
            
            q,gq = self.eval_Q(p,r,tol=tol,problem=problem)

            # Show progress
            if not quiet and i > 0:
                print '%d\t%.5e' %(i,Q)

            # Infinity
            if not q < np.inf:
                return np.inf,None

            # Update
            Q += (q-Q)/(i+1.)
            gQ += (gq-gQ)/(i+1.)
                            
        return Q,gQ

    def eval_EQ_parallel(self,p,samples=500,num_procs=None,tol=1e-4,quiet=True):
        """
        Evaluates E[Q(p,r)] and its gradient in parallel. 

        Parameters
        ----------
        p : generator powers
        samples : number of samples
        num_procs : number of parallel processes
        tol : evaluation tolerance
        quiet : flag
        """
    
        if not num_procs:
            num_procs = cpu_count()
        pool = Pool(num_procs)
        num = int(np.ceil(float(samples)/float(num_procs)))
        results = zip(*pool.map(ApplyFunc,[(self,'eval_EQ',p,num,i,tol,quiet) for i in range(num_procs)],chunksize=1))
        pool.terminate()
        pool.join()
        return map(lambda vals: sum(map(lambda val: val/float(num_procs),vals)),results)
        
    def eval_Q(self,p,r,quiet=True,check=False,tol=1e-4,problem=None,return_data=False):
        """
        Evaluates Q(p,r) and its gradient.

        Parameters
        ----------
        p : generator powers
        r : renewable powers
        quiet : flag
        check : flag
        tol : evaluation tolerance 
        problem : QuadProblem
        return_data : flag

        Returns
        -------
        Q : value
        gQ : gradient
        """

        # Local vars
        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r
        num_bus = self.num_bus
        num_br = self.num_br
        
        # Check
        assert(np.all(r > 0))
        
        # Problem
        if not problem:
            problem = self.get_problem_for_Q(p,r)
               
        try:

            # Solver
            solver = OptSolverIQP()
            solver.set_parameters({'quiet':quiet,
                                   'tol':tol})
            
            # Solve
            solver.solve(problem)

            # Info
            x = solver.get_primal_variables()
            lam,nu,mu,pi = solver.get_dual_variables()
            k = solver.get_iterations()
            q = x[:num_p]
            
            # Check
            if check:
                w = x[num_p:num_p+num_w]
                s = x[num_p+num_w:num_p+num_w+num_r]
                z = x[num_p+num_w+num_r:]
                assert(norm(self.G*(p+q)+self.R*s-self.A*w-self.b) < (1e-6)*norm(self.b))
                assert(norm(self.J*w-z) < (1e-6)*norm(z))
                assert(np.all(self.p_min <= p+q))
                assert(np.all(p+q <= self.p_max))
                assert(np.all(0 <= s))
                assert(np.all(s <= r))
                assert(np.all(self.z_min <= z))
                assert(np.all(self.z_max >= z))

            # Objective
            Q = 0.5*np.dot(q,self.H1*q)+np.dot(self.g1,q)
            
            # Gradient
            gQ = -(self.H1*q+self.g1)

            # Return
            if not return_data:
                return Q,gQ
            else:

                # Sensitivity (as linear operator)
                dqdpT = self.get_sol_sensitivity(problem,x,lam,mu,pi) # NEED TO GENERALIZE THIS TO dxdpT
                
                data = {'q': q,
                        'dqdpT': dqdpT}

                return Q,gQ,data

        # Errors
        except OptSolverError_MaxIters:
            return np.inf,None
        except OptSolverError_LineSearch:
            return np.inf,None
        except Exception,e:
            raise

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

    def get_strong_convexity_constant(self):
        """
        Gets strong convexity constant, which for
        this case is lambda_min(H0)/2.

        Returns
        -------
        c : float
        """

        H0 = self.H0.tocoo()
        assert(np.all(H0.row == H0.col))
        lmin = np.min(H0.data)
        assert(lmin > 0)
        return lmin/2.

    def get_problem_for_Q(self,p,r):
        """
        Constructs second-stage quadratic problem
        in the form

        minimize(x)   (1/2)x^THx + g^Tx
        subject to    Ax = b
                      l <= x <= u,

        where x = (q,w,s,z).

        Parameters
        ----------
        p : generator powers
        r : renewable powers

        Returns
        -------
        problem : QuadProblem
        """

        # Constatns
        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r
        num_bus = self.num_bus
        num_br = self.num_br
        Ow = coo_matrix((num_w,num_w))
        Os = coo_matrix((num_r,num_r))
        Oz = coo_matrix((num_br,num_br))
        Iz = eye(num_br,format='coo')
        ow = np.zeros(num_w)
        os = np.zeros(num_r)
        oz = np.zeros(num_br)
        cost_factor = self.parameters['cost_factor']

        H1 = self.H1/cost_factor
        g1 = self.g1/cost_factor
        
        # Form QP problem
        H = bmat([[H1,None,None,None],  # q: gen power adjustments
                  [None,Ow,None,None],  # w: bus voltage angles
                  [None,None,Os,None],  # s: curtailed renewable powers
                  [None,None,None,Oz]], # z: slack variables for thermal limits
                 format='coo')
        g = np.hstack((g1,ow,os,oz))
        A = bmat([[self.G,-self.A,self.R,None],
                  [None,self.J,None,-Iz]],format='coo')
        b = np.hstack((self.b-self.G*p,oz))
        l = np.hstack((self.p_min-p,
                       self.w_min,
                       os,
                       self.z_min))
        u = np.hstack((self.p_max-p,
                       self.w_max,
                       r,
                       self.z_max))
        
        # Return
        return QuadProblem(H,g,A,b,l,u)

    def eval_F(self,x,w,tol=1e-4):
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
        Q,gQ = self.eval_Q(x,w,tol=tol)

        return (phi+Q,gphi+gQ)

    def eval_F_approx(self,x,tol=1e-4):
        """
        Evaluates certainty equivalent
        version of objective function.

        Parameters
        ----------
        x : generator powers

        Returns
        -------
        F : float
        gF : vector
        """

        return self.eval_F(x,self.Er,tol=tol)

    def eval_EF(self,x,samples=500,tol=1e-4):
        """
        Evaluates expected objective function.

        Parameters
        ----------
        x : generator powers
        samples : number of samples
        
        Returns
        -------
        val : float
        """
        
        phi = 0.5*np.dot(x,self.H0*x)+np.dot(self.g0,x)
        gphi = self.H0*x + self.g0
        Q,gQ = self.eval_EQ_parallel(x,samples=samples,tol=tol)

        return (phi+Q,gphi+gQ)

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

    def save_x_info(self,x,filename):

        # Local variables
        num_p = self.num_p
        H0 = self.H0*np.ones(num_p)
        
        # Check
        assert(x.size == num_p)

        # Writer
        f = open(filename,'w')
        writer = csv.writer(f)
        
        # Write
        writer.writerow(['p','pmin','pmax','H0','g0'])
        for i in range(num_p):
            writer.writerow([x[i],
                             self.p_min[i],
                             self.p_max[i],
                             H0[i],
                             self.g0[i]])

        # Close
        f.close()

    def show(self):
        """
        Shows problem information.
        """

        Ctot = np.sum(self.r_max)
        Btot = np.sum(self.r_base)
        
        print 'Stochastic Two-Stage DCOPF'
        print '--------------------------'
        print 'buses            : %d' %self.num_bus
        print 'gens             : %d' %self.num_p
        print 'vargens          : %d' %self.num_r
        print 'penetration cap  : %.2f (%% of load)' %(100.*Ctot/self.total_load)
        print 'penetration base : %.2f (%% of load)' %(100.*Btot/self.total_load)
        print 'penetration std  : %.2f (%% of local cap)' %self.uncertainty
        print 'correlation rad  : %d (edges)' %(self.corr_radius)
        print 'correlation val  : %.2f (unitless)' %(self.corr_value)

    def solve_approx(self,g_corr=None,tol=1e-4,quiet=False,samples=500,init_data=None):
        """
        Solves certainty equivalent problem
        with slope correction.

        Parameters
        ----------
        g_corr : slope correction
        tol : optimality tolerance
        quiet : flag
        samples : number of samples

        Returns
        -------
        p : generator powers
        """
        
        # Constants
        num_p = self.num_p
        num_w = self.num_w
        num_r = self.num_r
        num_bus = self.num_bus
        num_br = self.num_br
        Ow = coo_matrix((num_w,num_w))
        Os = coo_matrix((num_r,num_r))
        Oz = coo_matrix((num_br,num_br))
        Onp = coo_matrix((num_bus,num_p))
        Onz = coo_matrix((num_bus,num_br))
        ow = np.zeros(num_w)
        os = np.zeros(num_r)
        oz = np.zeros(num_br)
        Iz = eye(num_br,format='coo')
        
        # Slope correction
        if g_corr is None:
            g_corr = np.zeros(num_p)

        # Problem construction
        H = bmat([[self.H0+self.H1,-self.H1,None,None,None], # p
                  [-self.H1,self.H1,None,None,None],         # y = p + q
                  [None,None,Ow,None,None],                  # w
                  [None,None,None,Os,None],                  # s
                  [None,None,None,None,Oz]],                 # z
                 format='coo')
        g = np.hstack((self.g0-self.g1+g_corr,self.g1,ow,os,oz))
        A = bmat([[Onp,self.G,-self.A,self.R,Onz],
                  [None,None,self.J,None,-Iz]],format='coo')
        b = np.hstack((self.b,oz))
        l = np.hstack((self.p_min,  # p
                       self.p_min,  # y
                       self.w_min,  # w
                       os,          # s
                       self.z_min)) # z
        u = np.hstack((self.p_max,
                       self.p_max,
                       self.w_max,
                       self.Er,
                       self.z_max))

        # Problem
        if init_data is None:
            problem = QuadProblem(H,g,A,b,l,u)
        else:
            problem = QuadProblem(H,g,A,b,l,u,
                                  x=init_data['x'],
                                  lam=init_data['lam'],
                                  mu=init_data['mu'],
                                  pi=init_data['pi'])
        
        # Solve
        solver = OptSolverIQP()
        solver.set_parameters({'quiet':quiet,
                               'tol':tol})
        solver.solve(problem)
        results = solver.get_results()
        x = results['x']
        lam = results['lam']
        mu = results['mu']
        pi = results['pi']

        # Check
        problem.eval(x)
        gphi = problem.gphi
        assert(norm(gphi-A.T*lam+mu-pi) < (1e-6)*(norm(gphi)+norm(lam)+norm(mu)+norm(pi)))
        assert(norm(mu*(u-x)) < (1e-6)*(norm(mu)+norm(u-x)))
        assert(norm(pi*(x-l)) < (1e-6)*(norm(pi)+norm(x-l)))
        assert(np.all(x < u))
        assert(np.all(x > l))
        assert(norm(A*x-b) < (1e-6)*norm(b))

        # Return
        return x[:num_p],results

    def get_sol_sensitivity(self,problem,x,lam,mu,pi):
        """
        Computes sensitivity of optimal q with respect to p
        as a linear operator.

        NEED TO EXTEND TO ALL x
        
        Parameters
        ----------
        problem : QuadProblem
        x : primal
        lam : dual var (eq constraints)
        mu : dual var (upper bounds)
        pi : dual var (lower bounds)
        
        Returns
        -------
        dqdpT : Linear Operator
        """

        H = problem.H
        g = problem.g
        A = problem.A
        b = problem.b
        u = problem.u
        l = problem.l

        n = A.shape[1]
        m = A.shape[0]

        Dmu = spdiags(mu,0,n,n)
        Dpi = spdiags(pi,0,n,n)
        Dux = spdiags(u-x,0,n,n)
        Dxl = spdiags(x-l,0,n,n)

        In = eye(n)
        
        K = bmat([[H,-A.T,In,-In],
                  [A,None,None,None],
                  [-Dmu,None,Dux,None],
                  [Dpi,None,None,Dxl]],
                 format='coo')
        KT = K.T
        
        Ibar = eye(self.num_p,K.shape[0])
        
        Onp = coo_matrix((n,self.num_p))
        bp = bmat([[-self.G],
                   [coo_matrix((self.num_br,self.num_p))]],
                  format='coo')
        up = -eye(n,self.num_p)
        lp = -eye(n,self.num_p)

        eta_p = bmat([[Onp],
                      [bp],
                      [-Dmu*up],
                      [Dpi*lp]],
                     format='coo')

        linsolver = new_linsolver('mumps','unsymmetric')
        linsolver.analyze(KT)
        linsolver.factorize(KT)
        
        dqdpT = LinearOperator((self.num_p,self.num_p),
                               lambda y : eta_p.T*linsolver.solve(Ibar.T*y))
        
        return dqdpT
