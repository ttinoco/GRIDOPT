#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import time
import numpy as np
from .method_error import *
from .method import PFmethod
from numpy.linalg import norm

class ACPF(PFmethod):
    """
    AC power flow method.
    """

    SHUNTS_MODE_LOCKED = 'locked'
    SHUNTS_MODE_FREE = 'free'
    SHUNTS_MODE_REG = 'regulating'

    name = 'ACPF'
    
    _parameters = {'weight_vmag': 1e0,         # weight for voltage magnitude regularization
                   'weight_vang': 1e-3,        # weight for angle difference regularization
                   'weight_powers': 1e-3,      # weight for gen powers regularization
                   'weight_controls': 1e0,     # weight for control deviation penalty
                   'weight_var': 1e-5,         # weight for generic regularization
                   'v_min_clip': 0.5,          # min v threshold for clipping
                   'v_max_clip': 1.5,          # max v threshold for clipping
                   'lin_pf': False,            # flag for using linearized power flow 
                   'limit_vars': True,         # flag for enforcing generator and VSC reactive power limits
                   'shunts_mode': 'locked',    # switched shunts mode: locked, free, regulating
                   'lock_taps': True,          # flag for locking transformer tap ratios
                   'lock_vsc_P_dc': True,      # flag for locking vsc P dc
                   'lock_csc_P_dc': True,      # flag for locking csc P dc
                   'lock_csc_i_dc': True,      # flag for locking csc i dc
                   'vdep_loads': False,        # flag for voltage dependent loads
                   'tap_step': 0.5,            # tap ratio acceleration factor (NR only)
                   'shunt_step': 0.5,          # susceptance acceleration factor (NR only)
                   'dtap': 1e-4,               # tap ratio perturbation (NR only)
                   'dsus': 1e-4,               # susceptance perturbation (NR only)
                   'pvpq_start_k': 0,          # start iteration number for PVPQ switching heuristics
                   'vmin_thresh': 0.1,         # threshold for vmin
                   'solver': 'nr'}             # OPTALG optimization solver (augl, ipopt, nr, inlp)

    _parameters_augl = {'feastol' : 1e-4,
                        'optol' : 1e0,
                        'kappa' : 1e-5,
                        'theta_max': 1e-6,
                        'sigma_init_max': 1e9}

    _parameters_ipopt = {}

    _parameters_inlp = {'feastol' : 1e-4,
                        'optol' : 1e0}
    
    _parameters_nr = {}
                  
    def __init__(self):

        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverNR, OptSolverINLP

        # Parent init
        PFmethod.__init__(self)

        # Solver params
        augl_params = OptSolverAugL.parameters.copy()
        augl_params.update(self._parameters_augl)   # overwrite defaults

        ipopt_params = OptSolverIpopt.parameters.copy()
        ipopt_params.update(self._parameters_ipopt) # overwrite defaults

        inlp_params = OptSolverINLP.parameters.copy()
        inlp_params.update(self._parameters_inlp)   # overwrite defaults

        nr_params = OptSolverNR.parameters.copy()
        nr_params.update(self._parameters_nr)       # overwrite defaults

        self._parameters = ACPF._parameters.copy()
        self._parameters['solver_parameters'] = {'augl': augl_params,
                                                 'ipopt': ipopt_params,
                                                 'nr': nr_params,
                                                 'inlp': inlp_params}

    def create_problem(self,net):

        import pfnet

        # Parameters
        params = self._parameters
        wm = params['weight_vmag']
        wa = params['weight_vang']
        wp = params['weight_powers']
        wc = params['weight_controls']
        wv = params['weight_var']
        lin_pf = params['lin_pf']
        limit_vars = params['limit_vars']
        shunts_mode = params['shunts_mode']
        lock_taps = params['lock_taps']
        lock_vsc_P_dc = params['lock_vsc_P_dc']
        lock_csc_P_dc = params['lock_csc_P_dc']
        lock_csc_i_dc = params['lock_csc_i_dc']
        vdep_loads = params['vdep_loads']
        solver_name = params['solver']

        # Shunts mode
        if shunts_mode not in [self.SHUNTS_MODE_LOCKED,
                               self.SHUNTS_MODE_FREE,
                               self.SHUNTS_MODE_REG]:
            raise ValueError('invalid shunts mode')

        # Linearized model
        if lin_pf:
            limit_vars = False
        
        # Clear flags
        net.clear_flags()

        # OPT-based
        ###########
        if solver_name != 'nr':

            # Buses
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')
            if not limit_vars:
                net.set_flags('bus',
                              'fixed',
                              'v set regulated',
                              'voltage magnitude')

            # Genertors
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')

            # Loads
            if vdep_loads:
                for load in net.loads:
                    if load.is_voltage_dependent() and load.is_in_service():
                        net.set_flags_of_component(load,
                                                   'variable',
                                                   ['active power', 'reactive power'])
            # VSC HVDC
            net.set_flags('vsc converter',
                          'variable',
                          'any',
                          ['dc power', 'active power', 'reactive power'])

            # CSC HVDC
            net.set_flags('csc converter',
                          'variable',
                          'any',
                          ['dc power', 'active power', 'reactive power'])

            # DC buses
            net.set_flags('dc bus',
                          'variable',
                          'any',
                          'voltage')

            # FACTS
            net.set_flags('facts',
                          'variable',
                          'any',
                          'all')

            # Tap changers
            if not lock_taps:
                net.set_flags('branch', 
                              'variable',
                              'tap changer - v',
                              'tap ratio')

            # Swtiched shunts
            if shunts_mode != self.SHUNTS_MODE_LOCKED:
                net.set_flags('shunt', 
                              'variable',
                              'switching - v',
                              'susceptance')

            # Checks
            try:
                num_vars = (2*net.get_num_buses(True)-net.get_num_slack_buses(True) +
                            net.get_num_dc_buses(True) +
                            net.get_num_slack_gens(True) +
                            net.get_num_reg_gens(True) +
                            net.get_num_vsc_converters(True)*4 +
                            net.get_num_csc_converters(True)*4 +
                            net.get_num_facts(True)*9)*net.num_periods
                if not lock_taps:
                    num_vars += net.get_num_tap_changers_v(True)*net.num_periods
                if shunts_mode != self.SHUNTS_MODE_LOCKED:
                    num_vars += net.get_num_switched_v_shunts(True)*net.num_periods
                if vdep_loads:
                    num_vars += 2*net.get_num_vdep_loads(True)*net.num_periods
                    
                assert(net.num_vars == num_vars)
                assert(net.num_bounded == 0)
                if limit_vars:
                    assert(net.num_fixed == 0)
                else:                    
                    assert(net.num_fixed == len([b for b in net.buses if b.is_v_set_regulated(True)])*net.num_periods)
            except AssertionError:
                raise PFmethodError_BadProblem()  
            
            # Set up problem
            problem = pfnet.Problem(net)

            if lin_pf:
                problem.add_constraint(pfnet.Constraint('linearized AC power balance', net))
            else:
                problem.add_constraint(pfnet.Constraint('AC power balance', net))
            problem.add_constraint(pfnet.Constraint('HVDC power balance', net))
            problem.add_constraint(pfnet.Constraint('generator active power participation', net))
            problem.add_constraint(pfnet.Constraint('VSC converter equations', net))
            problem.add_constraint(pfnet.Constraint('CSC converter equations', net))
            problem.add_constraint(pfnet.Constraint('VSC DC voltage control', net))
            problem.add_constraint(pfnet.Constraint('CSC DC voltage control', net))
            problem.add_constraint(pfnet.Constraint('power factor regulation', net))
            problem.add_constraint(pfnet.Constraint('FACTS equations', net))
            if lock_vsc_P_dc:
                problem.add_constraint(pfnet.Constraint('VSC DC power control', net))
            if lock_csc_P_dc:
                problem.add_constraint(pfnet.Constraint('CSC DC power control', net))
            if lock_csc_i_dc:
                problem.add_constraint(pfnet.Constraint('CSC DC current control', net))

            problem.add_function(pfnet.Function('variable regularization',
                                                wv/max([net.num_vars,1.]), net))
            problem.add_function(pfnet.Function('voltage magnitude regularization',
                                                wm/max([net.get_num_buses(True),1.]), net))
            problem.add_function(pfnet.Function('voltage angle regularization',
                                                wa/max([net.get_num_buses(True),1.]), net))
            problem.add_function(pfnet.Function('generator powers regularization',
                                                wp/max([net.get_num_generators(True),1.]), net))
            problem.add_function(pfnet.Function('VSC DC power control',
                                                wc/max([net.get_num_vsc_converters(True),1.]), net))
            problem.add_function(pfnet.Function('CSC DC power control',
                                                wc/max([net.get_num_csc_converters(True),1.]), net))
            problem.add_function(pfnet.Function('CSC DC current control',
                                                wc/max([net.get_num_csc_converters(True),1.]), net))
            problem.add_function(pfnet.Function('FACTS active power control',
                                                wc/max([net.num_facts,1.]), net))
            problem.add_function(pfnet.Function('FACTS reactive power control',
                                                wc/max([net.num_facts,1.]), net))
            if limit_vars:
                problem.add_constraint(pfnet.Constraint('voltage set point regulation', net))
            else:
                problem.add_constraint(pfnet.Constraint('variable fixing', net))
            if not lock_taps:
                problem.add_constraint(pfnet.Constraint('voltage regulation by transformers', net))
                problem.add_function(pfnet.Function('tap ratio regularization',
                                                    wc/max([net.get_num_tap_changers_v(True),1.]), net))
            if shunts_mode != self.SHUNTS_MODE_LOCKED:
                if shunts_mode == self.SHUNTS_MODE_REG:
                    problem.add_constraint(pfnet.Constraint('voltage regulation by shunts', net))
                problem.add_function(pfnet.Function('susceptance regularization',
                                                    wc/max([net.get_num_switched_v_shunts(True),1.]), net))
            if vdep_loads:
                problem.add_constraint(pfnet.Constraint('load voltage dependence', net))
            problem.analyze()

            # Return
            return problem

        # NR-based
        ##########
        elif solver_name == 'nr':

            # Buses
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          'voltage angle')
            net.set_flags('bus',
                          'variable',
                          'any',
                          'voltage magnitude')

            # DC buses
            net.set_flags('dc bus',
                          'variable',
                          'any',
                          'voltage')

            # Generators
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')

            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
            
            # VSC HVDC
            net.set_flags('vsc converter',
                          'variable',
                          'any',
                          ['dc power', 'active power', 'reactive power'])

            # CSC HVDC
            net.set_flags('csc converter',
                          'variable',
                          'any',
                          ['dc power', 'active power', 'reactive power'])
            
            # FACTS
            net.set_flags('facts',
                          'variable',
                          'any',
                          'all')
            
            # Loads
            if vdep_loads:
                for load in net.loads:
                    if load.is_voltage_dependent() and load.is_in_service():
                        net.set_flags_of_component(load,
                                                   'variable',
                                                   ['active power', 'reactive power'])

            # Tap changers
            net.set_flags('branch',
                          ['variable','fixed'],
                          'tap changer - v',
                          'tap ratio')

            # Switched shunts
            net.set_flags('shunt',
                          ['variable','fixed'],
                          'switching - v',
                          'susceptance')
                  
            # Checks
            try:
                num_vars = (2*net.get_num_buses(True)-net.get_num_slack_buses(True) +
                            net.get_num_dc_buses(True) +
                            net.get_num_slack_gens(True) +
                            net.get_num_reg_gens(True) +
                            net.get_num_tap_changers_v(True) +
                            net.get_num_switched_v_shunts(True) +
                            net.get_num_vsc_converters(True)*4 +
                            net.get_num_csc_converters(True)*4 +
                            net.get_num_facts(True)*9)*net.num_periods
                if vdep_loads:
                    num_vars += 2*net.get_num_vdep_loads(True)*net.num_periods
                    
                assert(net.num_vars == num_vars)
                assert(net.num_fixed == (net.get_num_tap_changers_v(True) +
                                         net.get_num_switched_v_shunts(True))*net.num_periods)
            except AssertionError:
                raise PFmethodError_BadProblem()

            # Set up problem
            problem = pfnet.Problem(net)

            if lin_pf:
                problem.add_constraint(pfnet.Constraint('linearized AC power balance', net))
            else:
                problem.add_constraint(pfnet.Constraint('AC power balance', net))
            problem.add_constraint(pfnet.Constraint('HVDC power balance', net))
            problem.add_constraint(pfnet.Constraint('generator active power participation', net))
            problem.add_constraint(pfnet.Constraint('variable fixing', net))
            problem.add_constraint(pfnet.Constraint('PVPQ switching', net))
            problem.add_constraint(pfnet.Constraint('VSC converter equations', net))
            problem.add_constraint(pfnet.Constraint('CSC converter equations', net))
            problem.add_constraint(pfnet.Constraint('VSC DC voltage control', net))
            problem.add_constraint(pfnet.Constraint('CSC DC voltage control', net))
            problem.add_constraint(pfnet.Constraint('VSC DC power control', net))
            problem.add_constraint(pfnet.Constraint('CSC DC power control', net))
            problem.add_constraint(pfnet.Constraint('CSC DC current control', net))
            problem.add_constraint(pfnet.Constraint('switching power factor regulation', net))
            problem.add_constraint(pfnet.Constraint('FACTS equations', net))
            problem.add_constraint(pfnet.Constraint('switching FACTS active power control', net))
            problem.add_constraint(pfnet.Constraint('switching FACTS reactive power control', net))
            if vdep_loads:
                problem.add_constraint(pfnet.Constraint('load voltage dependence', net))
            if limit_vars:
                problem.add_heuristic(pfnet.Heuristic('PVPQ switching', net))
                problem.add_heuristic(pfnet.Heuristic('switching power factor regulation', net))
            problem.analyze()

            # Check
            if (problem.J.shape[0] + problem.A.shape[0] !=
                problem.get_num_primal_variables()):
                raise PFmethodError_BadProblem()
            
            # Return
            return problem

        # Invalid
        #########
        else:
            raise PFmethodError_BadOptSolver()
            
    def solve(self, net, save_problem=False):

        from optalg.opt_solver import OptSolverError, OptTermination, OptCallback
        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverNR, OptSolverINLP
        
        # Parameters
        params = self._parameters
        shunts_mode = params['shunts_mode']
        lock_taps= params['lock_taps']
        vmin_thresh = params['vmin_thresh']
        solver_name = params['solver']
        solver_params = params['solver_parameters']
        v_min_clip = params['v_min_clip']
        v_max_clip = params['v_max_clip']

        # Opt solver
        if solver_name == 'augl':
            solver = OptSolverAugL()
        elif solver_name == 'ipopt':
            solver = OptSolverIpopt()
        elif solver_name == 'inlp':
            solver = OptSolverINLP()
        elif solver_name == 'nr':
            solver = OptSolverNR()
        else:
            raise PFmethodError_BadOptSolver()
        solver.set_parameters(solver_params[solver_name])

        # Copy network
        net = net.get_copy(merge_buses=True)
        self.set_network_snapshot(net)

        # Clipping
        for bus in net.buses:
            bus.v_mag = np.minimum(np.maximum(bus.v_mag, v_min_clip), v_max_clip)

        # Problem
        t0 = time.time()
        problem = self.create_problem(net)
        problem_time = time.time()-t0

        # Callbacks
        def c1(s):
            if (s.k != 0 and
                (not lock_taps) and norm(s.problem.f,np.inf) < 100.*solver_params['nr']['feastol']):
                try:
                    self.apply_tran_v_regulation(s)
                except Exception as e:
                    raise PFmethodError_TranVReg(e)
            
        def c2(s):
            if (s.k != 0 and
                (shunts_mode == self.SHUNTS_MODE_REG) and
                norm(s.problem.f,np.inf) < 100.*solver_params['nr']['feastol']):
                try:
                    self.apply_shunt_v_regulation(s)
                except Exception as e:
                    raise PFmethodError_ShuntVReg(e)                

        def c3(s):
            if s.k >= params['pvpq_start_k']:
                prob = s.problem.wrapped_problem
                prob.apply_heuristics(s.x)
                s.problem.A = prob.A
                s.problem.b = prob.b

        if solver_name == 'nr':
            solver.add_callback(OptCallback(c1))
            solver.add_callback(OptCallback(c2))
            solver.add_callback(OptCallback(c3))
                
        # Termination
        def t1(s):
            if np.min(s.problem.wrapped_problem.network.bus_v_min) < vmin_thresh:
                return True
            else:
                return False
        solver.add_termination(OptTermination(t1,'low voltage'))
            
        # Info printer
        info_printer = self.get_info_printer()
        solver.set_info_printer(info_printer)
        
        # Solve
        update = True
        t0 = time.time()
        try:
            solver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(e)
        except Exception as e:
            update = False
            raise e
        finally:
            
            # Update network
            if update:
                net.set_var_values(solver.get_primal_variables()[:net.num_vars])
                net.update_properties()
                net.clear_sensitivities()
                if solver_name != 'nr':
                    problem.store_sensitivities(*solver.get_dual_variables())

            # Save results
            self.set_solver_name(solver_name)
            self.set_solver_status(solver.get_status())
            self.set_solver_message(solver.get_error_msg())
            self.set_solver_iterations(solver.get_iterations())
            self.set_solver_time(time.time()-t0)
            self.set_solver_primal_variables(solver.get_primal_variables())
            self.set_solver_dual_variables(solver.get_dual_variables())
            self.set_problem(problem if save_problem else None)
            self.set_problem_time(problem_time)
            self.set_network_snapshot(net)
 
    def get_info_printer(self):

        # Parameters
        solver_name = self._parameters['solver']

        # OPT-based
        ###########
        if solver_name != 'nr':
        
            def info_printer(solver,header):
                net = solver.problem.wrapped_problem.network
                if header:
                    print('{0:^5}'.format('vmax'), end=' ')
                    print('{0:^5}'.format('vmin'), end=' ')
                    print('{0:^8}'.format('gvdev'), end=' ')
                    print('{0:^8}'.format('gQvio'))
                else:
                    print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                    print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                    print('{0:^8.1e}'.format(np.average(net.gen_v_dev)), end=' ')
                    print('{0:^8.1e}'.format(np.average(net.gen_Q_vio)))
            return info_printer

        # NR-based
        ##########
        elif solver_name == 'nr':

            def info_printer(solver,header):
                net = solver.problem.wrapped_problem.network
                if header:
                    print('{0:^5}'.format('vmax'), end=' ')
                    print('{0:^5}'.format('vmin'), end=' ')
                    print('{0:^8}'.format('gvdev'), end=' ')
                    print('{0:^8}'.format('gQvio'))
                else:
                    print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                    print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                    print('{0:^8.1e}'.format(np.average(net.gen_v_dev)), end=' ')
                    print('{0:^8.1e}'.format(np.average(net.gen_Q_vio)))
            return info_printer

        # Invalid
        #########
        else:
            raise PFmethodError_BadOptSolver()

    def apply_shunt_v_regulation(self,solver):

        # Local variables
        dsus = self._parameters['dsus']
        step = self._parameters['shunt_step']
        p = solver.problem.wrapped_problem
        net = p.network
        x = solver.x
        eps = 1e-8

        # Fix constraints
        c = p.find_constraint('variable fixing')
        A = c.A        
        b = c.b
        
        # Rhs
        rhs = np.hstack((np.zeros(p.f.size),np.zeros(p.b.size)))

        # Offset
        offset = 0
        for c in p.constraints:
            if c.name == 'variable fixing':
                break
            else:
                offset += c.A.shape[0]

        # Violation check
        for i in range(net.num_buses):

            bus = net.get_bus(i)

            if bus.is_regulated_by_shunt(True) and not bus.is_slack():
                
                assert(bus.has_flags('variable','voltage magnitude'))

                for t in range(net.num_periods):

                    v = x[bus.index_v_mag[t]]
                    vmax = bus.v_max_reg
                    vmin = bus.v_min_reg

                    assert(len(bus.reg_shunts) > 0)
                    assert(vmax >= vmin)

                    # Violation
                    if v > vmax or v < vmin:

                        for reg_shunt in bus.reg_shunts:

                            if not reg_shunt.is_in_service():
                                continue
                    
                            assert(reg_shunt.has_flags('variable','susceptance'))
                        
                            s = x[reg_shunt.index_b[t]]
                            smax = reg_shunt.b_max
                            smin = reg_shunt.b_min
                            assert(smin <= smax)
                            
                            # Fix constr index
                            k = int(np.where(A.col == reg_shunt.index_b[t])[0])
                            i = A.row[k]
                            
                            assert(np.abs(b[i]-x[reg_shunt.index_b[t]]) < eps)
                            assert(A.data[k] == 1.)
                            
                            # Sensitivity
                            assert(rhs[p.f.size+offset+i] == 0.)
                            rhs[p.f.size+offset+i] = dsus
                            dx = solver.linsolver.solve(rhs)
                            dv = dx[bus.index_v_mag[t]]
                            dvds = dv/dsus
                            rhs[p.f.size+offset+i] = 0.
                            
                            # Adjustment
                            dv = (vmax+vmin)/2.-v
                            ds = step*dv/dvds if dvds != 0. else 0.
                            snew = np.maximum(np.minimum(s+ds,smax),smin)
                            x[reg_shunt.index_b[t]] = snew
                            b[i] = snew
                            if np.abs(snew-s) > eps:
                                break
                            
        # Update
        solver.func(x)
        p.update_lin()
        solver.problem.A = p.A
        solver.problem.b = p.b

    def apply_tran_v_regulation(self,solver):
        
        # Local variables
        dtap = self._parameters['dtap']
        step = self._parameters['tap_step']
        p = solver.problem.wrapped_problem
        net = p.network
        x = solver.x
        eps = 1e-8

        # Fix constraints
        c = p.find_constraint('variable fixing')
        A = c.A        
        b = c.b

        # Rhs
        rhs = np.hstack((np.zeros(p.f.size),np.zeros(p.b.size)))
        
        # Offset
        offset = 0
        for c in p.constraints:
            if c.name == 'variable fixing':
                break
            else:
                offset += c.A.shape[0]

        # Violation check
        for i in range(net.num_buses):

            bus = net.get_bus(i)

            if bus.is_regulated_by_tran(True) and not bus.is_slack():
                
                assert(bus.has_flags('variable','voltage magnitude'))

                for tau in range(net.num_periods):
                    
                    v = x[bus.index_v_mag[tau]]
                    vmax = bus.v_max_reg
                    vmin = bus.v_min_reg
                    
                    assert(len(bus.reg_trans) > 0)
                    assert(vmax > vmin)
                
                    # Violation
                    if v > vmax or v < vmin:
                        
                        for reg_tran in bus.reg_trans:

                            if not reg_tran.is_in_service():
                                continue
                                                        
                            assert(reg_tran.has_flags('variable','tap ratio'))
                            
                            t = x[reg_tran.index_ratio[tau]]
                            tmax = reg_tran.ratio_max
                            tmin = reg_tran.ratio_min
                            assert(tmax >= tmin)
                            
                            # Fix constr index
                            k = int(np.where(A.col == reg_tran.index_ratio[tau])[0])
                            i = A.row[k]
                            assert(np.abs(b[i]-x[reg_tran.index_ratio[tau]]) < eps)
                            assert(A.data[k] == 1.)
                            
                            # Sensitivity
                            assert(rhs[p.f.size+offset+i] == 0.)
                            rhs[p.f.size+offset+i] = dtap
                            dx = solver.linsolver.solve(rhs)
                            dv = dx[bus.index_v_mag[tau]]
                            dvdt = dv/dtap
                            rhs[p.f.size+offset+i] = 0.
                            
                            # Adjustment
                            dv = (vmax+vmin)/2.-v
                            dt = step*dv/dvdt if dvdt != 0. else 0.
                            tnew = np.maximum(np.minimum(t+dt,tmax),tmin)
                            x[reg_tran.index_ratio[tau]] = tnew
                            b[i] = tnew
                            if np.abs(tnew-t) > eps:
                                break

        # Update
        solver.func(x)        
        p.update_lin()
        solver.problem.A = p.A
        solver.problem.b = p.b
