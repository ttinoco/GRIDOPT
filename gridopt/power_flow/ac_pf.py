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
    
    CONTROL_MODE_LOCKED = 'locked'
    CONTROL_MODE_FREE = 'free'
    CONTROL_MODE_REG = 'regulating'

    name = 'ACPF'
    
    _parameters = {'weight_vmag': 1e0,       # weight for voltage magnitude regularization
                   'weight_vang': 1e-3,      # weight for angle difference regularization
                   'weight_powers': 1e-3,    # weight for gen powers regularization
                   'weight_controls': 1e0,   # weight for control deviation penalty
                   'weight_var': 1e-5,       # weight for general variable regularization
                   'v_min_clip': 0.5,        # lower v threshold for clipping
                   'v_max_clip': 1.5,        # upper v threshold for clipping
                   'v_limits': False,        # voltage magnitude limits
                   'Q_limits': True,         # flag for enforcing generator, VSC and FACTS reactive power limits
                   'Q_mode': 'regulating',   # reactive power mode: free, regulating
                   'shunt_limits': True,     # flag for enforcing switched shunt susceptance limits
                   'shunt_mode': 'locked',   # switched shunts mode: locked, free, regulating
                   'tap_limits': True,       # flag for enforcing transformer tap ratio limits
                   'tap_mode': 'locked',     # transformer tap ratio mode: locked, free, regulating
                   'lock_vsc_P_dc': True,    # flag for locking vsc P dc
                   'lock_csc_P_dc': True,    # flag for locking csc P dc
                   'lock_csc_i_dc': True,    # flag for locking csc i dc
                   'vdep_loads': False,      # flag for modeling voltage dependent loads
                   'pvpq_start_k': 0,        # start iteration number for PVPQ switching heuristics
                   'vmin_thresh': 0.1,       # minimum voltage magnitude threshold
                   'gens_redispatch': False, # flag for allowing active power redispatch
                   'shunts_round': True,     # flag for rounding discrete switched shunt susceptances (not supported yet)
                   'taps_round': True,       # flag for rounding discrete transformer tap ratios (not supported yet)
                   'v_mag_warm_ref': False,  # flag for using current v_mag as reference in v_mag regularization
                   'solver': 'nr',           # OPTALG optimization solver: augl, ipopt, nr, inlp
                   'tap_step': 0.5,          # tap ratio acceleration factor (NR only)
                   'shunt_step': 0.5,        # susceptance acceleration factor (NR only)
                   'dtap': 1e-4,             # tap ratio perturbation (NR only)
                   'dsus': 1e-4}             # susceptance perturbation (NR only)

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

    def create_problem(self, net):

        solver_name = self._parameters['solver']

        if solver_name == 'nr':
            return self.create_problem_nr(net)
        else:
            return self.create_problem_opt(net)

    def create_problem_nr(self, net):

        import pfnet

        # Parameters
        params = self._parameters
        Q_mode = params['Q_mode']
        Q_limits = params['Q_limits']
        shunt_mode = params['shunt_mode']
        shunt_limits = params['shunt_limits']
        tap_mode = params['tap_mode']
        tap_limits = params['tap_limits']
        lock_vsc_P_dc = params['lock_vsc_P_dc']
        lock_csc_P_dc = params['lock_csc_P_dc']
        lock_csc_i_dc = params['lock_csc_i_dc']
        vdep_loads = params['vdep_loads']
        gens_redispatch = params['gens_redispatch']

        # Check shunt options 
        if shunt_mode not in [self.CONTROL_MODE_LOCKED,
                               self.CONTROL_MODE_REG]:
            raise ValueError('invalid shunts mode')
        if shunt_mode == self.CONTROL_MODE_REG and not shunt_limits:
            raise ValueError('unsupported shunts configuration')

        # Check tap options 
        if tap_mode not in [self.CONTROL_MODE_LOCKED,
                               self.CONTROL_MODE_REG]:
            raise ValueError('invalid taps mode')
        if tap_mode == self.CONTROL_MODE_REG and not tap_limits:
            raise ValueError('unsupported taps configuration')
        
        # Check Q options
        if Q_mode != self.CONTROL_MODE_REG:
            raise ValueError('invalid reactive power mode')
        
        # Check other options
        if gens_redispatch:
            raise ValueError('generation redispatch not supported')
        if not lock_vsc_P_dc:
            raise ValueError('VSC P DC must be locked')
        if not lock_csc_P_dc:
            raise ValueError('CSC P DC must be locked')
        if not lock_csc_i_dc:
            raise ValueError('CSC i DC must be locked')
        
        # Clear flags
        net.clear_flags()

        # Buses
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')
        net.set_flags('bus',
                      'variable',
                      'any',
                      'voltage magnitude')

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
        
        # Loads
        if vdep_loads:
            for load in net.loads:
                if load.is_voltage_dependent() and load.is_in_service():
                    net.set_flags_of_component(load,
                                               'variable',
                                               ['active power', 'reactive power'])
                    
        # Tap changers 
        if tap_mode != self.CONTROL_MODE_LOCKED:
            net.set_flags('branch',
                          ['variable', 'fixed'],
                          'tap changer - v',
                          'tap ratio')
        
        # Switched shunts
        if shunt_mode != self.CONTROL_MODE_LOCKED:
            net.set_flags('shunt',
                          ['variable', 'fixed'],
                          'switching - v',
                          'susceptance')
        
        # Set up problem
        problem = pfnet.Problem(net)
        
        problem.add_constraint(pfnet.Constraint('AC power balance', net))
        problem.add_constraint(pfnet.Constraint('HVDC power balance', net))
        problem.add_constraint(pfnet.Constraint('generator active power participation', net))
        problem.add_constraint(pfnet.Constraint('variable fixing', net))
        problem.add_constraint(pfnet.Constraint('VSC converter equations', net))
        problem.add_constraint(pfnet.Constraint('CSC converter equations', net))
        problem.add_constraint(pfnet.Constraint('FACTS equations', net))
        problem.add_constraint(pfnet.Constraint('VSC DC voltage control', net))
        problem.add_constraint(pfnet.Constraint('CSC DC voltage control', net))
        problem.add_constraint(pfnet.Constraint('VSC DC power control', net))
        problem.add_constraint(pfnet.Constraint('CSC DC power control', net))
        problem.add_constraint(pfnet.Constraint('CSC DC current control', net))

        problem.add_constraint(pfnet.Constraint('PVPQ switching', net))
        problem.add_constraint(pfnet.Constraint('switching power factor regulation', net))
        problem.add_constraint(pfnet.Constraint('switching FACTS active power control', net))
        problem.add_constraint(pfnet.Constraint('switching FACTS reactive power control', net))
        
        if vdep_loads:
            problem.add_constraint(pfnet.Constraint('load voltage dependence', net))

        if Q_limits:            
            problem.add_heuristic(pfnet.Heuristic('PVPQ switching', net))
            problem.add_heuristic(pfnet.Heuristic('switching power factor regulation', net))

        problem.analyze()
        
        # Check
        if (problem.J.shape[0] + problem.A.shape[0] != problem.get_num_primal_variables()):
            raise PFmethodError_BadProblem()
        
        # Return
        return problem

    def create_problem_opt(self, net):

        import pfnet

        # Parameters
        params = self._parameters
        wm = params['weight_vmag']
        wa = params['weight_vang']
        wp = params['weight_powers']
        wc = params['weight_controls']
        wv = params['weight_var']
        v_limits = params['v_limits']
        Q_mode = params['Q_mode']
        Q_limits = params['Q_limits']
        shunt_mode = params['shunt_mode']
        shunt_limits = params['shunt_limits']
        tap_mode = params['tap_mode']
        tap_limits = params['tap_limits']
        lock_vsc_P_dc = params['lock_vsc_P_dc']
        lock_csc_P_dc = params['lock_csc_P_dc']
        lock_csc_i_dc = params['lock_csc_i_dc']
        vdep_loads = params['vdep_loads']
        v_mag_warm_ref = params['v_mag_warm_ref']
        gens_redispatch = params['gens_redispatch']
        
        # Check shunt options
        if shunt_mode not in [self.CONTROL_MODE_LOCKED,
                               self.CONTROL_MODE_FREE,
                               self.CONTROL_MODE_REG]:
            raise ValueError('invalid shunts mode')
        if shunt_mode == self.CONTROL_MODE_REG and not shunt_limits:
            raise ValueError('unsupported shunts configuration')

        # Check tap options
        if tap_mode not in [self.CONTROL_MODE_LOCKED,
                             self.CONTROL_MODE_FREE,
                             self.CONTROL_MODE_REG]:
            raise ValueError('invalid taps mode')
        if tap_mode == self.CONTROL_MODE_REG and not tap_limits:
            raise ValueError('unsupported taps configuration')

        # Check Q options
        if Q_mode not in [self.CONTROL_MODE_REG,
                          self.CONTROL_MODE_FREE]:
            raise ValueError('invalid reactive power mode')

        # Clear flags
        net.clear_flags()

        # Buses
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')
        net.set_flags('bus',
                      'variable',
                      'any',
                      'voltage magnitude')
        if Q_mode == self.CONTROL_MODE_REG and not Q_limits: 
            net.set_flags('bus',
                          'fixed',
                          'v set regulated',
                          'voltage magnitude')
        if v_limits:
            net.set_flags('bus',
                          'bounded',
                          'any',
                          'voltage magnitude')
            
        # Genertors 
        if gens_redispatch:    
            net.set_flags('generator',
                          ['variable', 'bounded'],
                          'any',
                          'active power')
        else:  
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
        net.set_flags('generator',
                      'variable',
                      'regulator',
                      'reactive power')
        if Q_mode == self.CONTROL_MODE_FREE and Q_limits: 
            net.set_flags('generator',
                          'bounded',
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
        if Q_mode == self.CONTROL_MODE_FREE and Q_limits: 
            net.set_flags('vsc converter',
                          'bounded',
                          'any',
                          'reactive power')
            
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
        if Q_mode == self.CONTROL_MODE_FREE and Q_limits:        
            net.set_flags('facts',
                          'bounded',
                          'any',
                          'reactive power')
            
        # Tap changers 
        if tap_mode != self.CONTROL_MODE_LOCKED:
            net.set_flags('branch', 
                          'variable',
                          'tap changer - v',
                          'tap ratio')
        if tap_mode == self.CONTROL_MODE_FREE and tap_limits:
            net.set_flags('branch',
                          'bounded',
                          'tap changer - v',
                          'tap ratio')
            
        # Swtiched shunts
        if shunt_mode != self.CONTROL_MODE_LOCKED:
            net.set_flags('shunt', 
                          'variable',
                          'switching - v',
                          'susceptance')
        if shunt_mode == self.CONTROL_MODE_FREE and shunt_limits:
            net.set_flags('shunt',
                          'bounded',
                          'switching - v',
                          'susceptance')
            
        # Set up problem
        problem = pfnet.Problem(net)
        
        problem.add_constraint(pfnet.Constraint('AC power balance', net))
        problem.add_constraint(pfnet.Constraint('HVDC power balance', net))
        problem.add_constraint(pfnet.Constraint('generator active power participation', net))
        problem.add_constraint(pfnet.Constraint('VSC converter equations', net))
        problem.add_constraint(pfnet.Constraint('CSC converter equations', net))
        problem.add_constraint(pfnet.Constraint('FACTS equations', net))
        problem.add_constraint(pfnet.Constraint('VSC DC voltage control', net))
        problem.add_constraint(pfnet.Constraint('CSC DC voltage control', net))
        problem.add_constraint(pfnet.Constraint('power factor regulation', net))
        
        if lock_vsc_P_dc:
            problem.add_constraint(pfnet.Constraint('VSC DC power control', net))
        if lock_csc_P_dc:
            problem.add_constraint(pfnet.Constraint('CSC DC power control', net))
        if lock_csc_i_dc:
            problem.add_constraint(pfnet.Constraint('CSC DC current control', net))

        func = pfnet.Function('voltage magnitude regularization', wm/(net.get_num_buses(True)+1.), net)
        func.set_parameter('v_set_reference', not v_mag_warm_ref)
        problem.add_function(func)

        problem.add_function(pfnet.Function('variable regularization', wv/(net.num_vars+1.), net))
        problem.add_function(pfnet.Function('voltage angle regularization', wa/(net.get_num_buses(True)+1.), net))
        problem.add_function(pfnet.Function('generator powers regularization', wp/(net.get_num_generators(True)+1.), net))
        problem.add_function(pfnet.Function('VSC DC power control', wc/(net.get_num_vsc_converters(True)+1.), net))
        problem.add_function(pfnet.Function('CSC DC power control', wc/(net.get_num_csc_converters(True)+1.), net))
        problem.add_function(pfnet.Function('CSC DC current control', wc/(net.get_num_csc_converters(True)+1.), net))
        problem.add_function(pfnet.Function('FACTS active power control', wc/(net.get_num_facts(True)+1.), net))
        problem.add_function(pfnet.Function('FACTS reactive power control', wc/(net.get_num_facts(True)+1.), net))
        
        if gens_redispatch:
            problem.add_function(pfnet.Function('generation redispatch penalty', wc/(net.get_num_generators(True)+1.), net))
            
        if Q_mode == self.CONTROL_MODE_REG and Q_limits: 
            problem.add_constraint(pfnet.Constraint('voltage set point regulation', net))

        if net.num_fixed > 0:
            problem.add_constraint(pfnet.Constraint('variable fixing', net))
            
        if tap_mode != self.CONTROL_MODE_LOCKED:
            problem.add_function(pfnet.Function('tap ratio regularization', wc/(net.get_num_tap_changers_v(True)+1.), net))
            if tap_mode == self.CONTROL_MODE_REG and tap_limits: 
                problem.add_constraint(pfnet.Constraint('voltage regulation by transformers', net))
            
        if shunt_mode != self.CONTROL_MODE_LOCKED:
            problem.add_function(pfnet.Function('susceptance regularization', wc/(net.get_num_switched_v_shunts(True)+1.), net))
            if shunt_mode == self.CONTROL_MODE_REG and shunt_limits:
                problem.add_constraint(pfnet.Constraint('voltage regulation by shunts', net))
            
        if vdep_loads:
            problem.add_constraint(pfnet.Constraint('load voltage dependence', net))
            
        if net.num_bounded > 0: 
            problem.add_constraint(pfnet.Constraint('variable bounds', net))
            
        # Analyze
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self, net, save_problem=False):

        from optalg.opt_solver import OptSolverError, OptTermination, OptCallback
        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverNR, OptSolverINLP
        
        # Parameters
        params = self._parameters
        Q_mode = params['Q_mode']
        shunt_mode = params['shunt_mode']
        shunts_round = params['shunts_round']
        tap_mode = params['tap_mode']
        taps_round = params['taps_round']
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
            if (s.k != 0 and params['tap_limits'] and tap_mode == self.CONTROL_MODE_REG and
                norm(s.problem.f, np.inf) < 100.*solver_params['nr']['feastol']):
                try:
                    self.apply_tran_v_regulation(s)
                except Exception as e:
                    raise PFmethodError_TranVReg(e)
            
        def c2(s):
            if (s.k != 0 and params['shunt_limits'] and shunt_mode == self.CONTROL_MODE_REG and
                norm(s.problem.f, np.inf) < 100.*solver_params['nr']['feastol']):
                try:
                    self.apply_shunt_v_regulation(s)
                except Exception as e:
                    raise PFmethodError_ShuntVReg(e)                

        def c3(s):
            if (s.k >= params['pvpq_start_k'] and params['Q_limits'] and Q_mode == self.CONTROL_MODE_REG):
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
        solver.add_termination(OptTermination(t1, 'low voltage'))
            
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

        # Define
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

        # Return
        return info_printer

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

