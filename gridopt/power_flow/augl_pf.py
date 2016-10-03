#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import pfnet
import numpy as np
from .method_error import *
from .method import PFmethod
from optalg.opt_solver import OptSolverError,OptCallback,OptTermination,OptSolverAugL

class AugLPF(PFmethod):
    """
    Augmented Lagrangian-based power flow method.
    """

    name = 'AugLPF'
        
    parameters = {'weight_vmag':1e0,     # for reg voltage magnitude penalty
                  'weight_vang':1e-3,    # for angle difference penalty
                  'weight_pq':1e-3,      # for gen powers
                  'weight_t':1e1,        # for tap ratios
                  'weight_b':1e-4,       # for shunt susceptances
                  'lock_taps':True,      # flag for locking transformer tap ratios
                  'lock_shunts':True,    # flag for locking swtiched shunts
                  'vmin_thresh':0.1}     # threshold for vmin
                  
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = AugLPF.parameters.copy()
        self.parameters.update(OptSolverAugL.parameters)

    def create_problem(self,net):

        # Parameters
        params = self.parameters
        weight_vang = params['weight_vang']
        weight_vmag = params['weight_vmag']
        weight_pq = params['weight_pq']
        weight_t = params['weight_t']
        weight_b = params['weight_b']
        lock_taps = params['lock_taps']
        lock_shunts = params['lock_shunts']
        
        # Clear flags
        net.clear_flags()

        # Adjust gens
        net.adjust_generators()
        
        # Set up variables
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VMAG|pfnet.BUS_VAR_VANG)
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK|pfnet.BUS_PROP_REG_BY_GEN,
                      pfnet.BUS_VAR_VDEV)
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_SLACK,
                      pfnet.GEN_VAR_P)
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS,
                      pfnet.GEN_PROP_REG,
                      pfnet.GEN_VAR_Q)
        
        # Tap ratios
        if not lock_taps:
            net.set_flags(pfnet.OBJ_BUS,
                          pfnet.FLAG_VARS,
                          pfnet.BUS_PROP_REG_BY_TRAN,
                          pfnet.BUS_VAR_VVIO)
            net.set_flags(pfnet.OBJ_BRANCH, 
                          pfnet.FLAG_VARS,
                          pfnet.BRANCH_PROP_TAP_CHANGER_V,
                          pfnet.BRANCH_VAR_RATIO|pfnet.BRANCH_VAR_RATIO_DEV)

        # Shunt voltage control
        if not lock_shunts:
            net.set_flags(pfnet.OBJ_BUS,
                          pfnet.FLAG_VARS,
                          pfnet.BUS_PROP_REG_BY_SHUNT,
                          pfnet.BUS_VAR_VVIO)
            net.set_flags(pfnet.OBJ_SHUNT, 
                          pfnet.FLAG_VARS,
                          pfnet.SHUNT_PROP_SWITCHED_V,
                          pfnet.SHUNT_VAR_SUSC|pfnet.SHUNT_VAR_SUSC_DEV)

        try:
            num_vars = (2*(net.num_buses-net.get_num_slack_buses()) +
                        2*(net.get_num_buses_reg_by_gen()-net.get_num_slack_buses()) +
                        net.get_num_slack_gens() +
                        net.get_num_reg_gens())*net.num_periods
            if not lock_taps:
                num_vars += (2*net.get_num_buses_reg_by_tran() +
                             3*net.get_num_tap_changers_v())*net.num_periods
            if not lock_shunts:
                num_vars += (2*net.get_num_buses_reg_by_shunt() +
                             3*net.get_num_switched_shunts())*net.num_periods
            if not lock_shunts and not lock_taps:
                num_vars -= 2*len([b for b in net.buses 
                                   if (b.is_regulated_by_tran() and
                                       b.is_regulated_by_shunt())])*net.num_periods
            assert(net.num_vars == num_vars)
        except AssertionError:
            raise PFmethodError_BadProblem(self)  
        
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_PF)
        problem.add_constraint(pfnet.CONSTR_TYPE_PAR_GEN_P)
        problem.add_constraint(pfnet.CONSTR_TYPE_PAR_GEN_Q)
        problem.add_constraint(pfnet.CONSTR_TYPE_REG_GEN)
        if not lock_taps:
            problem.add_constraint(pfnet.CONSTR_TYPE_REG_TRAN)
        if not lock_shunts:
            problem.add_constraint(pfnet.CONSTR_TYPE_REG_SHUNT)
        problem.add_function(pfnet.FUNC_TYPE_REG_VMAG,weight_vmag)
        problem.add_function(pfnet.FUNC_TYPE_REG_VANG,weight_vang)
        problem.add_function(pfnet.FUNC_TYPE_REG_PQ,weight_pq)
        if not lock_taps:
            problem.add_function(pfnet.FUNC_TYPE_REG_RATIO,weight_t)
        if not lock_shunts:
            problem.add_function(pfnet.FUNC_TYPE_REG_SUSC,weight_b)
        problem.analyze()
        
        # Return
        return problem

    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.network
            if header:
                print('{0:^5}'.format('vmax'), end=' ')
                print('{0:^5}'.format('vmin'), end=' ')
                print('{0:^8}'.format('gvdev'), end=' ')
                print('{0:^8}'.format('gQvio'), end=' ')
                print('{0:^8}'.format('tvvio'), end=' ')
                print('{0:^8}'.format('svvio'))
            else:
                print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.gen_v_dev)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.gen_Q_vio)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.tran_v_vio)), end=' ')
                print('{0:^8.1e}'.format(np.average(net.shunt_v_vio)))
        return info_printer
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        vmin_thresh = params['vmin_thresh']

        # Problem
        problem = self.create_problem(net)

        # Callbacks

        # Termination
        def t1(s):
            if np.min(s.problem.network.bus_v_min) < vmin_thresh:
                return True
            else:
                return False
        
        # Info printer
        info_printer = self.get_info_printer()
            
        # Set up solver
        solver = OptSolverAugL()
        solver.set_parameters(params)
        solver.add_termination(OptTermination(t1,'low voltage'))
        solver.set_info_printer(info_printer)
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError as e:
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
        assert(problem.A.shape[0] == lam.size)
        assert(problem.f.shape[0] == nu.size)
        assert(mu is None or not mu.size)
        assert(pi is None or not pi.size)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam,nu,mu,pi)
