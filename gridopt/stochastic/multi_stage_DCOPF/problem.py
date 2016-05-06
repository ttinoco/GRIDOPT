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
from optalg.opt_solver.opt_solver_error import *
from optalg.stoch_solver import StochObjMS_Problem
from optalg.opt_solver import OptSolverIQP,QuadProblem
from scipy.sparse import triu,bmat,coo_matrix,eye,spdiags
            
class MS_DCOPF(StochObjMS_Problem):
    
    # Parameters
    parameters = {'cost_factor' : 1e2,   # factor for determining fast gen cost
                  'infinity' : 1e3,      # infinity
                  'flow_factor' : 1.,    # factor for relaxing thermal limits
                  'max_ramping' : 0.1,   # factor for constructing ramping limits
                  'num_samples' : 2000}  # number of samples

    def __init__(self,net,horizon):
        """
        Class constructor.
        
        Parameters
        ----------
        net : PFNET Network
        """

        # Parameters
        self.parameters = MS_DCOPF.parameters.copy()

        # Save info
        self.net = net
        self.T = horizon

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
        num_l = net.num_loads                           # loads
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
        net.set_flags(pf.OBJ_LOAD,
                      pf.FLAG_VARS,
                      pf.LOAD_PROP_ANY,
                      pf.LOAD_VAR_P)
        net.set_flags(pf.OBJ_VARGEN,
                      pf.FLAG_VARS,
                      pf.VARGEN_PROP_ANY,
                      pf.VARGEN_VAR_P)
        assert(net.num_vars == num_w+num_p+num_r+num_l)

        # Current values
        x = net.get_var_values()

        # Projections
        Pw = net.get_var_projection(pf.OBJ_BUS,pf.BUS_VAR_VANG)
        Pp = net.get_var_projection(pf.OBJ_GEN,pf.GEN_VAR_P)
        Pl = net.get_var_projection(pf.OBJ_LOAD,pf.LOAD_VAR_P)
        Pr = net.get_var_projection(pf.OBJ_VARGEN,pf.VARGEN_VAR_P)
        assert(Pw.shape == (num_w,net.num_vars))
        assert(Pp.shape == (num_p,net.num_vars))
        assert(Pl.shape == (num_l,net.num_vars))
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
        H = (cost.Hphi + cost.Hphi.T - triu(cost.Hphi))/net.base_power # symmetric, scaled
        g = cost.gphi/net.base_power - H*x
        
        # Bounds
        l = net.get_var_values(pf.LOWER_LIMITS)
        u = net.get_var_values(pf.UPPER_LIMITS)
        assert(np.all(Pw*l < Pw*u))
        assert(np.all(Pp*l < Pp*u))
        assert(np.all(Pl*l <= Pl*u))
        assert(np.all(Pr*l < Pr*u))

        # Problem data
        self.num_p = num_p
        self.num_q = num_p
        self.num_w = num_w
        self.num_s = num_r
        self.num_r = num_r
        self.num_y = num_p
        self.num_z = num_br
        self.num_l = num_l
        self.num_bus = num_bus
        self.num_br = num_br

        self.p_max = Pp*u
        self.p_min = Pp*l
        
        self.q_max = Pp*u
        self.q_min = Pp*l
        
        self.w_max = self.parameters['infinity']*np.ones(self.num_w)
        self.w_min = -self.parameters['infinity']*np.ones(self.num_w)

        self.r_max = Pr*u
        self.r_base = Pr*x

        self.z_max = hu
        self.z_min = hl

        self.y_max = self.parameters['max_ramping']*(self.p_max-self.p_min)
        self.y_min = -self.parameters['max_ramping']*(self.p_max-self.p_min)
        
    def show(self):
        """
        Shows problem information.
        """

        tot_cap = np.sum(self.r_max)
        tot_base = np.sum(self.r_base)
        tot_unc = sum([g.P_std for g in self.net.var_generators])
        tot_load = sum([l.P for l in self.net.loads])
        
        print 'Stochastic Multi-Stage DCOPF'
        print '----------------------------'
        print 'buses            : %d' %self.num_bus
        print 'gens             : %d' %self.num_p
        print 'vargens          : %d' %self.num_r
        print 'time horizon     : %d' %self.T
        print 'penetration cap  : %.2f (%% of load)' %(100.*tot_cap/tot_load)
        print 'penetration base : %.2f (%% of load)' %(100.*tot_base/tot_load)
        print 'penetration std  : %.2f (%% of local cap)' %(100.*tot_unc/tot_load)
        print 'correlation rad  : %d (edges)' %(self.net.vargen_corr_radius)
        print 'correlation val  : %.2f (unitless)' %(self.net.vargen_corr_value)        
