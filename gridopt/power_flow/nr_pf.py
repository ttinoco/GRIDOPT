#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from method_error import *
from method import PFmethod
from optalg.opt_solver import OptSolverError,OptCallback,OptTermination,OptSolverNR

class NRPF(PFmethod):
    """
    Newton-Raphson power flow method.
    """

    name = 'NRPF'
        
    parameters = {'limit_gens':True,  # enforce generator reactive power limits
                  'lock_taps':True,   # lock tap ratios of transformer
                  'tap_step':0.5,     # tap ratio acceleration factor
                  'lock_shunts':True, # lock susceptances of shunt devices
                  'shunt_step':0.5,   # susceptance acceleration factor
                  'vmin_thresh':0.1,  # low voltage threshold
                  'dtap':1e-5,        # tap ratio perturbation
                  'dsus':1e-5}        # susceptance perturbation
    
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = NRPF.parameters.copy()       # method parameters
        self.parameters.update(OptSolverNR.parameters) # solver parameters
        
    def apply_shunt_v_regulation(self,solver):

        # Local variables
        dsus = self.parameters['dsus']
        step = self.parameters['shunt_step']
        p = solver.problem
        net = p.network
        x = solver.x
        eps = 1e-8

        # Fix constraints
        c = p.find_constraint(pfnet.CONSTR_TYPE_FIX)
        A = c.A        
        b = c.b

        # Rhs
        rhs = np.hstack((np.zeros(p.f.size),np.zeros(p.b.size)))

        # Offset
        offset = 0
        for c in p.constraints:
            if c.type == pfnet.CONSTR_TYPE_FIX:
                break
            else:
                offset += c.A.shape[0]

        # Violation check
        for i in range(net.num_buses):

            bus = net.get_bus(i)

            if bus.is_regulated_by_shunt() and not bus.is_slack():

                assert(bus.has_flags(pfnet.FLAG_VARS,pfnet.BUS_VAR_VMAG))

                v = x[bus.index_v_mag]
                vmax = bus.v_max
                vmin = bus.v_min

                assert(len(bus.reg_shunts) > 0)
                assert(vmax >= vmin)

                # Violation
                if v > vmax or v < vmin:

                    for reg_shunt in bus.reg_shunts:
                    
                        assert(reg_shunt.has_flags(pfnet.FLAG_VARS,pfnet.SHUNT_VAR_SUSC))
                                
                        s = x[reg_shunt.index_b]
                        smax = reg_shunt.b_max
                        smin = reg_shunt.b_min
                        assert(smin <= smax)
                        
                        # Fix constr index
                        k = int(np.where(A.col == reg_shunt.index_b)[0])
                        i = A.row[k]
                        assert(np.abs(b[i]-x[reg_shunt.index_b]) < eps)
                        assert(A.data[k] == 1.)
                        
                        # Sensitivity
                        assert(rhs[p.f.size+offset+i] == 0.)
                        rhs[p.f.size+offset+i] = dsus
                        dx = solver.linsolver.solve(rhs)
                        dv = dx[bus.index_v_mag]
                        dvds = dv/dsus
                        rhs[p.f.size+offset+i] = 0.
                        
                        # Adjustment
                        dv = (vmax+vmin)/2.-v
                        ds = step*dv/dvds if dvds != 0. else 0.
                        snew = np.maximum(np.minimum(s+ds,smax),smin)
                        x[reg_shunt.index_b] = snew
                        b[i] = snew
                        if np.abs(snew-s) > eps:
                            break

        # Update
        solver.func(x)        
        solver.problem.update_lin()

    def apply_tran_v_regulation(self,solver):
        
        # Local variables
        dtap = self.parameters['dtap']
        step = self.parameters['tap_step']
        p = solver.problem
        net = p.network
        x = solver.x
        eps = 1e-8
            
        # Fix constraints
        c = p.find_constraint(pfnet.CONSTR_TYPE_FIX)
        A = c.A        
        b = c.b

        # Rhs
        rhs = np.hstack((np.zeros(p.f.size),np.zeros(p.b.size)))
        
        # Offset
        offset = 0
        for c in p.constraints:
            if c.type == pfnet.CONSTR_TYPE_FIX:
                break
            else:
                offset += c.A.shape[0]

        # Violation check
        for i in range(net.num_buses):

            bus = net.get_bus(i)

            if bus.is_regulated_by_tran() and not bus.is_slack():

                assert(bus.has_flags(pfnet.FLAG_VARS,pfnet.BUS_VAR_VMAG))

                v = x[bus.index_v_mag]
                vmax = bus.v_max
                vmin = bus.v_min

                assert(len(bus.reg_trans) > 0)
                assert(vmax > vmin)
                
                # Violation
                if v > vmax or v < vmin:

                    for reg_tran in bus.reg_trans:
                    
                        assert(reg_tran.has_flags(pfnet.FLAG_VARS,pfnet.BRANCH_VAR_RATIO))
                                
                        t = x[reg_tran.index_ratio]
                        tmax = reg_tran.ratio_max
                        tmin = reg_tran.ratio_min
                        assert(tmax >= tmin)
                        
                        # Fix constr index
                        k = int(np.where(A.col == reg_tran.index_ratio)[0])
                        i = A.row[k]
                        assert(np.abs(b[i]-x[reg_tran.index_ratio]) < eps)
                        assert(A.data[k] == 1.)
                    
                        # Sensitivity
                        assert(rhs[p.f.size+offset+i] == 0.)
                        rhs[p.f.size+offset+i] = dtap
                        dx = solver.linsolver.solve(rhs)
                        dv = dx[bus.index_v_mag]
                        dvdt = dv/dtap
                        rhs[p.f.size+offset+i] = 0.
                        
                        # Adjustment
                        dv = (vmax+vmin)/2.-v
                        dt = step*dv/dvdt if dvdt != 0. else 0.
                        tnew = np.maximum(np.minimum(t+dt,tmax),tmin)
                        x[reg_tran.index_ratio] = tnew
                        b[i] = tnew
                        if np.abs(tnew-t) > eps:
                            break

        # Update
        solver.func(x)        
        solver.problem.update_lin()

    def create_problem(self,net):
        
        # Parameters
        params = self.parameters
        limit_gens = params['limit_gens']
        
        # Clear flags
        net.clear_flags()

        # Adjust gens
        net.adjust_generators()
        
        # Voltages
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VMAG|pfnet.BUS_VAR_VANG)
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_FIXED,
                      pfnet.BUS_PROP_REG_BY_GEN,
                      pfnet.BUS_VAR_VMAG)

        # Gen active powers
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_SLACK,
                      pfnet.GEN_VAR_P)

        # Gen reactive powers
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_REG,
                      pfnet.GEN_VAR_Q)

        # Tap ratios
        net.set_flags(pfnet.OBJ_BRANCH,
                      pfnet.FLAG_VARS|pfnet.FLAG_FIXED,
                      pfnet.BRANCH_PROP_TAP_CHANGER_V,
                      pfnet.BRANCH_VAR_RATIO)

        # Shunt susceptances
        net.set_flags(pfnet.OBJ_SHUNT,
                      pfnet.FLAG_VARS|pfnet.FLAG_FIXED,
                      pfnet.SHUNT_PROP_SWITCHED_V,
                      pfnet.SHUNT_VAR_SUSC)

        try:
            assert(net.num_vars == (2*(net.num_buses-net.get_num_slack_buses()) +
                                    net.get_num_slack_gens() +
                                    net.get_num_reg_gens() +
                                    net.get_num_tap_changers() +
                                    net.get_num_switched_shunts()))
            assert(net.num_fixed == (net.get_num_buses_reg_by_gen() +
                                     net.get_num_tap_changers() +
                                     net.get_num_switched_shunts()))
        except AssertionError:
            raise PFmethodError_BadProblem(self)

        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_PF)
        problem.add_constraint(pfnet.CONSTR_TYPE_PAR_GEN_P)
        problem.add_constraint(pfnet.CONSTR_TYPE_PAR_GEN_Q)
        problem.add_constraint(pfnet.CONSTR_TYPE_FIX)
        if limit_gens:
            problem.add_heuristic(pfnet.HEUR_TYPE_PVPQ)
        problem.analyze()

        # Return
        return problem
        
    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.network
            if header:
                print '{0:^5}'.format('vmax'),
                print '{0:^5}'.format('vmin'),
                print '{0:^8}'.format('gvdev'),
                print '{0:^8}'.format('gQvio'),
                print '{0:^8}'.format('tvvio'),
                print '{0:^8}'.format('svvio')
            else:
                print '{0:^5.2f}'.format(net.bus_v_max),
                print '{0:^5.2f}'.format(net.bus_v_min),
                print '{0:^8.1e}'.format(net.gen_v_dev),
                print '{0:^8.1e}'.format(net.gen_Q_vio),
                print '{0:^8.1e}'.format(net.tran_v_vio),
                print '{0:^8.1e}'.format(net.shunt_v_vio)
        return info_printer
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        feastol = params['feastol']
        lock_taps= params['lock_taps']
        lock_shunts = params['lock_shunts']
        vmin_thresh = params['vmin_thresh']
        
        # Problem
        problem = self.create_problem(net)        

        # Callbacks
        def c1(s):
            if (s.k != 0 and
                (not lock_taps) and 
                np.linalg.norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_tran_v_regulation(s)
            
        def c2(s):
            if (s.k != 0 and
                (not lock_shunts) and 
                np.linalg.norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_shunt_v_regulation(s)

        def c3(s):
            if s.k > 0:
                s.problem.apply_heuristics(s.x)

        # Termination
        def t1(s):
            if s.problem.network.bus_v_min < vmin_thresh:
                return True
            else:
                return False

        # Info printer
        info_printer = self.get_info_printer()

        # Set up solver
        solver = OptSolverNR()
        solver.set_parameters(params)
        solver.add_callback(OptCallback(c1))
        solver.add_callback(OptCallback(c2))
        solver.add_callback(OptCallback(c3))
        solver.add_termination(OptTermination(t1,'low voltage'))
        solver.set_info_printer(info_printer)
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError,e:
            raise PFmethodError_SolverError(self,e)
        finally:

            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(solver.get_primal_variables())
            self.set_dual_variables(solver.get_dual_variables())
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)
    
    def update_network(self,net):
        
        # Get data
        problem = self.results['problem']
        x = self.results['primal_variables']
        lam,nu,mu,pi = self.results['dual_variables']
       
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == x.shape)
        assert(net.num_vars == x.size)
        assert(lam is None)
        assert(nu is None)
        assert(mu is None)
        assert(pi is None)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()

