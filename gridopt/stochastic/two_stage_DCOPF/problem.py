#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet as pf
import numpy as np
from numpy.linalg import norm
from utils import ApplyFunc, Cache
from optalg.opt_solver import OptSolverIQP, QuadProblem
from optalg.lin_solver import new_linsolver
from optalg.opt_solver.opt_solver_error import *
from multiprocessing import Pool, cpu_count
from scipy.sparse import triu,bmat,coo_matrix,eye,spdiags
from optalg.stoch_solver import StochProblem
            
class TwoStageDCOPF(StochProblem):
    """"
    This class represents a problem of the form
    
    minimize(p)   \varphi_0(p) + E_r[Q(p,r)]
    subject to    p_min <= p <= p_max
    
    where Q(p,r) is the optimal value of

    minimize(q,w,s,z)  \varphi_1(q)
    subjcet to         G(p+q) + Rs - Aw = b
                       Jw - z = 0
                       p_min <= p+q <= p_max
                       z_min <= z <= z_max
                       0 <= s <= r
    """

    # Problem constants
    P_MIN = 1e-5
    COST_FACTOR = 100.
    VANG_LIMIT = 1000.
    FLOW_LIMIT = 1000.
    FLOW_FACTOR = 1.

    def __init__(self,net,penetration=100.,uncertainty=25.,corr_radius=1,corr_value=0.1):
        """
        Class constructor.
        
        Parameters
        ----------
        net : PFNET Network
        penetration : float (% of load)
        uncertainty : float (% of load)
        corr_radius : float (number of edges)
        corr_value : float (correlation coefficient)
        """

        # Configuration
        self.penetration = penetration
        self.uncertainty = uncertainty
        self.corr_radius = corr_radius
        self.corr_value = corr_value
        assert(0. <= penetration)
        assert(0. <= uncertainty)
        assert(0 <= corr_radius and isinstance(corr_radius,int))
        assert(-1. <= corr_value <= 1.)
        
        # Renewable sources
        assert(net.num_vargens == 0)
        gen_buses = net.get_gen_buses()
        self.total_load = sum([l.P for l in net.loads])
        vargen_array = pf.VarGeneratorArray(len(gen_buses))
        net.set_vargen_array(vargen_array)
        net.set_vargen_buses(gen_buses)
        assert(net.num_vargens == len([b for b in net.buses if b.gens]))
        for vg in net.var_generators:
            vg.P_max = self.total_load/net.num_vargens
            vg.P_std = (uncertainty/100.)*vg.P_max
            assert(vg.P_max > 0)
            assert(vg.P_std > 0)
        assert(np.abs(sum([vg.P_max for vg in net.var_generators])-self.total_load) < 1e-10)
        
        # Generator limits
        for gen in net.generators:
            gen.P_min = 0.
            gen.P_max = np.maximum(gen.P_max,0.)
            assert(gen.P_min <= gen.P_max)

        # Branch flow limits
        for br in net.branches:
            if br.ratingA == 0.:
                br.ratingA = self.FLOW_LIMIT
            else:
                br.ratingA *= self.FLOW_FACTOR
        
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
        hl = fl_lim.hl.copy()
        hu = fl_lim.hu.copy()
        assert(np.all(hl < hu))
        
        # Generation cost
        cost = pf.Function(pf.FUNC_TYPE_GEN_COST,1.,net)
        cost.analyze()
        cost.eval(x)
        H = cost.Hphi + cost.Hphi.T - triu(cost.Hphi) # symmetric
        g = cost.gphi - H*x
        
        # Bounds
        l = net.get_var_values(pf.LOWER_LIMITS)
        u = net.get_var_values(pf.UPPER_LIMITS)
        assert(np.all(Pw*l < Pw*u))
        assert(np.all(Pp*l < Pp*u))
        assert(np.all(Pr*l < Pr*u))
        
        # Problem data
        self.num_w = num_w
        self.num_p = num_p
        self.num_r = num_r
        self.num_bus = num_bus
        self.num_br = num_br
        self.p_max = Pp*u
        self.p_min = Pp*l
        self.w_max = self.VANG_LIMIT*np.ones(self.num_w)
        self.w_min = -self.VANG_LIMIT*np.ones(self.num_w)
        self.r_max = Pr*u
        self.z_max = hu
        self.z_min = hl 
        self.H0 = Pp*H*Pp.T
        self.g0 = Pp*g
        self.H1 = self.H0*self.COST_FACTOR
        self.g1 = np.zeros(num_p)
        self.G = A*Pp.T
        self.R = A*Pr.T
        self.A = -A*Pw.T
        self.J = G*Pw.T
        self.b = b
        self.Pp = Pp
        self.Pw = Pw
        self.Pr = Pr

        # Rrenewable base
        assert(0 <= penetration)
        fraction = penetration*np.sum([l.P for l in net.loads])/(100.*np.sum(self.r_max))
        self.r_base = fraction*self.r_max

        # Renewable covariance
        from scikits.sparse.cholmod import cholesky
        r_cov = Pr*net.create_vargen_P_sigma(corr_radius,corr_value)*Pr.T
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
        
        # Checks
        assert(np.all(self.p_min == 0))
        assert(np.all(0 < self.r_max))
        assert(np.all(self.z_min < self.z_max))
        assert(np.all(self.p_min < self.p_max))
        assert(np.all(cost.Hphi.row == cost.Hphi.col))
        assert(np.all(cost.Hphi.data > 0))
        assert(norm(self.A.T*np.ones(self.num_bus)) < (1e-10)*np.sqrt(self.num_bus*1.))
        
    def eval_EQ(self,p,feastol=1e-4,samples=500,quiet=True):
        """
        Evaluates E[Q(p,r)].
        """
        
        # Local vars
        Q = 0.
        gQ = np.zeros(self.num_p)
        
        # Problem
        problem = self.get_problem_for_Q(p,np.ones(self.num_r))
        
        # Sampling loop
        np.random.seed()
        for i in range(samples):
            
            r = self.sample_w()

            problem.u[self.num_p+self.num_w:self.num_p+self.num_w+self.num_r] = r # Important
            
            q,gq = self.eval_Q(p,r,feastol=feastol,problem=problem)

            if not quiet and i > 0:
                print '%d\t%.5e' %(i,Q)

            if not q < np.inf:
                return np.inf, None

            Q += (q-Q)/(i+1.)
            gQ += (gq-gQ)/(i+1.)
                            
        return Q,gQ

    def eval_EQ_parallel(self,p,feastol=1e-4,samples=500,quiet=True,num_procs=None):
    
        if not num_procs:
            num_procs = cpu_count()
        pool = Pool(num_procs)
        num = int(np.ceil(float(samples)/float(num_procs)))
        results = zip(*pool.map(ApplyFunc,num_procs*[(self,'eval_EQ',p,feastol,num,quiet)]))
        return map(lambda vals: sum(map(lambda val: num*val/float(num*num_procs),vals)), results)
        
    def eval_Q(self,p,r,quiet=True,check=False,feastol=1e-4,problem=None):
        """
        Evaluates Q(p,r).

        Parameters
        ----------
        p : generator powers
        r : renewable powers
        quiet : flag
        check : flag
        problem : QuadProblem

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
                                   'feastol':feastol})
            
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
                assert(norm(self.G*(p+q)+self.R*s-self.A*w-self.b) < (1e-3)*norm(self.b))
                assert(norm(self.J*w-z) < (1e-3)*norm(z))
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
            return Q,gQ

        # Errors
        except OptSolverError_MaxIters:
            return np.inf, None
        except OptSolverError_LineSearch:
            return np.inf, None
        except Exception,e:
            raise

    def get_Ew(self,samples=500):

        r = 0
        for i in range(samples):
            r += (self.sample_w()-r)/(i+1)
        return r

    def get_size_x(self):

        return self.num_p

    def get_strong_convexity_constant(self):

        H0 = self.H0.tocoo()
        assert(np.all(H0.row == H0.col))
        lmin = np.min(H0.data)
        assert(lmin > 0)
        return lmin/2.

    def get_problem_for_Q(self,p,r):

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

        H1 = self.H1/self.COST_FACTOR
        g1 = self.g1
        
        # Form QP problem
        H = bmat([[H1,None,None,None],
                  [None,Ow,None,None],
                  [None,None,Os,None],
                  [None,None,None,Oz]],
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

    def eval_F(self,x,w):
        
        phi = 0.5*np.dot(x,self.H0*x)+np.dot(self.g0,x)
        gphi = self.H0*x + self.g0
        Q,gQ = self.eval_Q(x,w)

        return (phi+Q,gphi+gQ)

    def eval_EF(self,x,samples=500):
        
        phi = 0.5*np.dot(x,self.H0*x)+np.dot(self.g0,x)
        gphi = self.H0*x + self.g0
        Q,gQ = self.eval_EQ_parallel(x,samples=samples)

        return (phi+Q,gphi+gQ)

    def project_on_X(self,x):

        return np.maximum(np.minimum(x,self.p_max),self.p_min)

    def sample_w(self):

        return np.minimum(np.maximum(self.r_base+self.L*np.random.randn(self.num_r),1e-3),self.r_max)

    def show(self):

        Ctot = np.sum(self.r_max)
        Rtot = np.sum(self.r_base)
        
        print 'Stochastic Two-Stage DCOPF'
        print '--------------------------'
        print 'buses           : %d' %self.num_bus
        print 'gens            : %d' %self.num_p
        print 'vargens         : %d' %self.num_r
        print 'penetration max : %.2f (%% of load)' %(100.*Ctot/self.total_load)
        print 'penetration ave : %.2f (%% of load)' %(100.*Rtot/self.total_load)
        print 'penetration std : %.2f (%% of local cap)' %(self.uncertainty)
        print 'correlation rad : %d (edges)' %(self.corr_radius)
        print 'correlation val : %.2f (unitless)' %(self.corr_value)

    def solve_certainty_equivalent(self,g_corr=None,Ew=None,feastol=1e-4,quiet=False,samples=500):
        
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

        # Sample mean
        if Ew is None:
            Ew = self.get_Ew(samples=samples)

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
                       Ew,
                       self.z_max))

        # Problem
        problem = QuadProblem(H,g,A,b,l,u)
        
        # Solver
        solver = OptSolverIQP()
        solver.set_parameters({'quiet':quiet,
                               'feastol':feastol})

        # Solve
        solver.solve(problem)
        
        x = solver.get_primal_variables()
        p = x[:num_p]
        y = x[num_p:2*num_p]
        w = x[2*num_p:2*num_p+num_w]
        s = x[2*num_p+num_w:2*num_p+num_w+num_r]
        z = x[2*num_p+num_w+num_r:]
        
        # Check
        assert(np.all(self.z_min < self.J*w))
        assert(np.all(self.z_max > self.J*w))
        assert(norm(self.J*w-z) < 1e-6)

        # Return
        return x[:num_p]
