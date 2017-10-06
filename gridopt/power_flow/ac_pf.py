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

    name = 'ACPF'
    
    _parameters = {'weight_vmag': 1e0,  # weight for reg voltage magnitude penalty
                   'weight_vang': 1e0,  # weight for angle difference penalty
                   'weight_pq': 1e-3,   # weight for gen powers penalty
                   'weight_t': 1e-3,    # weight for tap ratios penalty
                   'weight_b': 1e-3,    # weight for shunt susceptances penalty
                   'limit_gens': True,  # flag for enforcing generator reactive power limits
                   'lock_taps': True,   # flag for locking transformer tap ratios
                   'lock_shunts': True, # flag for locking swtiched shunts
                   'tap_step': 0.5,     # tap ratio acceleration factor (NR only)
                   'shunt_step': 0.5,   # susceptance acceleration factor (NR only)
                   'dtap': 1e-5,        # tap ratio perturbation (NR only)
                   'dsus': 1e-5,        # susceptance perturbation (NR only)
                   'vmin_thresh': 0.1,  # threshold for vmin
                   'optsolver': 'augl'} # OPTALG optimization solver (augl,ipopt,nr)

    _parameters_augl = {'feastol' : 1e-4,
                        'optol' : 1e-4,
                        'kappa' : 1e-5}

    _parameters_ipopt = {}

    _parameters_inlp = {}
    
    _parameters_nr = {}
                  
    def __init__(self):

        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverNR, OptSolverINLP

        # Parent init
        PFmethod.__init__(self)

        # Optsolver params
        augl_params = OptSolverAugL.parameters.copy()
        augl_params.update(self._parameters_augl)   # overwrite defaults

        ipopt_params = OptSolverIpopt.parameters.copy()
        ipopt_params.update(self._parameters_ipopt) # overwrite defaults

        inlp_params = OptSolverINLP.parameters.copy()
        inlp_params.update(self._parameters_inlp)   # overwrite defaults

        nr_params = OptSolverNR.parameters.copy()
        nr_params.update(self._parameters_nr)       # overwrite defaults

        self._parameters.update(ACPF._parameters)
        self._parameters['optsolver_parameters'] = {'augl': augl_params,
                                                    'ipopt': ipopt_params,
                                                    'nr': nr_params,
                                                    'inlp': inlp_params}

    def create_problem(self,net):

        import pfnet

        # Parameters
        params = self._parameters
        wm = params['weight_vang']
        wa = params['weight_vmag']
        wp = params['weight_pq']
        wt = params['weight_t']
        wb = params['weight_b']
        limit_gens = params['limit_gens']
        lock_taps = params['lock_taps']
        lock_shunts = params['lock_shunts']
        optsolver_name = params['optsolver']
        
        # Clear flags
        net.clear_flags()

        # OPT-based
        ###########
        if optsolver_name != 'nr':
            
            # Set up variables
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage angle','voltage magnitude'])
            if not limit_gens:
                net.set_flags('bus',
                              'fixed',
                              'regulated by generator',
                              'voltage magnitude')
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')
        
            # Tap ratios
            if not lock_taps:
                net.set_flags('branch', 
                              'variable',
                              'tap changer - v',
                              'tap ratio')

            # Shunt voltage control
            if not lock_shunts:
                net.set_flags('shunt', 
                              'variable',
                              'switching - v',
                              'susceptance')

            try:
                num_vars = (2*(net.num_buses-net.get_num_slack_buses()) +
                            net.get_num_slack_gens() +
                            net.get_num_reg_gens())*net.num_periods
                if not lock_taps:
                    num_vars += net.get_num_tap_changers_v()*net.num_periods
                if not lock_shunts:
                    num_vars += net.get_num_switched_shunts()*net.num_periods
                assert(net.num_vars == num_vars)
                if limit_gens:
                    assert(net.num_fixed == 0)
                else:
                    assert(net.num_fixed == net.get_num_buses_reg_by_gen()*net.num_periods)
            except AssertionError:
                raise PFmethodError_BadProblem()  
            
            # Set up problem
            problem = pfnet.Problem(net)
            problem.add_constraint(pfnet.Constraint('AC power balance',net))
            problem.add_constraint(pfnet.Constraint('generator active power participation',net))
            problem.add_constraint(pfnet.Constraint('generator reactive power participation',net))
            problem.add_function(pfnet.Function('voltage magnitude regularization',
                                                wm/max([net.num_buses,1.]),net))
            problem.add_function(pfnet.Function('voltage angle regularization',
                                                wa/max([net.num_buses,1.]),net))
            problem.add_function(pfnet.Function('generator powers regularization',
                                                wp/max([net.num_generators,1.]),net))
            if limit_gens:
                problem.add_constraint(pfnet.Constraint('voltage regulation by generators',net))
            else:
                problem.add_constraint(pfnet.Constraint('variable fixing',net))
            if not lock_taps:
                problem.add_constraint(pfnet.Constraint('voltage regulation by transformers',net))
                problem.add_function(pfnet.Function('tap ratio regularization',
                                                    wt/max([net.get_num_tap_changers_v(),1.]),net))
            if not lock_shunts:
                problem.add_constraint(pfnet.Constraint('voltage regulation by shunts',net))
                problem.add_function(pfnet.Function('susceptance regularization',
                                                    wb/max([net.get_num_switched_shunts(),1.]),net))
            problem.analyze()
        
            # Return
            return problem

        # NR-based
        ##########
        elif optsolver_name == 'nr':

            # Voltages
            net.set_flags('bus',
                          'variable',
                          'not slack',
                          ['voltage magnitude','voltage angle'])
            net.set_flags('bus',
                          'fixed',
                          'regulated by generator',
                          'voltage magnitude')

            # Gen active powers
            net.set_flags('generator',
                          'variable',
                          'slack',
                          'active power')

            # Gen reactive powers
            net.set_flags('generator',
                          'variable',
                          'regulator',
                          'reactive power')

            # Tap ratios
            net.set_flags('branch',
                          ['variable','fixed'],
                          'tap changer - v',
                          'tap ratio')

            # Shunt susceptances
            net.set_flags('shunt',
                          ['variable','fixed'],
                          'switching - v',
                          'susceptance')

            try:
                assert(net.num_vars == (2*(net.num_buses-net.get_num_slack_buses()) +
                                        net.get_num_slack_gens() +
                                        net.get_num_reg_gens() +
                                        net.get_num_tap_changers_v() +
                                        net.get_num_switched_shunts())*net.num_periods)
                assert(net.num_fixed == (net.get_num_buses_reg_by_gen() +
                                         net.get_num_tap_changers_v() +
                                         net.get_num_switched_shunts())*net.num_periods)
            except AssertionError:
                raise PFmethodError_BadProblem()

            # Set up problem
            problem = pfnet.Problem(net)
            problem.add_constraint(pfnet.Constraint('AC power balance',net))
            problem.add_constraint(pfnet.Constraint('generator active power participation',net))
            problem.add_constraint(pfnet.Constraint('generator reactive power participation',net))
            problem.add_constraint(pfnet.Constraint('variable fixing',net))
            if limit_gens:
                problem.add_heuristic(pfnet.HEUR_TYPE_PVPQ)
            problem.analyze()

            # Return
            return problem

        # Invalid
        #########
        else:
            raise PFmethodError_BadOptSolver()
            
    def solve(self,net):

        from optalg.opt_solver import OptSolverError, OptTermination, OptCallback
        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverNR, OptSolverINLP
        
        # Parameters
        params = self._parameters
        lock_taps= params['lock_taps']
        lock_shunts = params['lock_shunts']
        vmin_thresh = params['vmin_thresh']
        optsolver_name = params['optsolver']
        optsolver_params = params['optsolver_parameters']
        feastol = optsolver_params['nr']['feastol']

        # Opt solver
        if optsolver_name == 'augl':
            optsolver = OptSolverAugL()
        elif optsolver_name == 'ipopt':
            optsolver = OptSolverIpopt()
        elif optsolver_name == 'inlp':
            optsolver = OptSolverINLP()
        elif optsolver_name == 'nr':
            optsolver = OptSolverNR()
        else:
            raise PFmethodError_BadOptSolver()
        optsolver.set_parameters(optsolver_params[optsolver_name])

        # Copy network
        net = net.get_copy()

        # Problem
        t0 = time.time()
        problem = self.create_problem(net)
        problem_time = time.time()-t0
        
        # Callbacks
        def c1(s):
            if (s.k != 0 and
                (not lock_taps) and 
                norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_tran_v_regulation(s)
            
        def c2(s):
            if (s.k != 0 and
                (not lock_shunts) and 
                norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_shunt_v_regulation(s)

        def c3(s):
            if s.k > 0:
                prob = s.problem.wrapped_problem
                prob.apply_heuristics(s.x)
                s.problem.A = prob.A
                s.problem.b = prob.b

        if optsolver_name == 'nr':
            optsolver.add_callback(OptCallback(c1))
            optsolver.add_callback(OptCallback(c2))
            optsolver.add_callback(OptCallback(c3))
                
        # Termination
        def t1(s):
            if np.min(s.problem.wrapped_problem.network.bus_v_min) < vmin_thresh:
                return True
            else:
                return False
        optsolver.add_termination(OptTermination(t1,'low voltage'))
            
        # Info printer
        info_printer = self.get_info_printer()
        optsolver.set_info_printer(info_printer)
        
        # Solve
        update = True
        t0 = time.time()
        try:
            optsolver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(e)
        except Exception as e:
            update = False
            raise e
        finally:
            
            # Update network
            if update:
                net.set_var_values(optsolver.get_primal_variables()[:net.num_vars])
                net.update_properties()
                net.clear_sensitivities()
                if optsolver_name != 'nr':
                    problem.store_sensitivities(*optsolver.get_dual_variables())

            # Save results
            self.set_solver_name(optsolver_name)
            self.set_solver_status(optsolver.get_status())
            self.set_solver_message(optsolver.get_error_msg())
            self.set_solver_iterations(optsolver.get_iterations())
            self.set_solver_time(time.time()-t0)
            self.set_solver_primal_variables(optsolver.get_primal_variables())
            self.set_solver_dual_variables(optsolver.get_dual_variables())
            self.set_problem(None) # skip for now
            self.set_problem_time(problem_time)
            self.set_network_snapshot(net)
 
    def get_info_printer(self):

        # Parameters
        optsolver_name = self._parameters['optsolver']

        # OPT-based
        ###########
        if optsolver_name != 'nr':
        
            def info_printer(solver,header):
                net = solver.problem.wrapped_problem.network
                if header:
                    print('{0:^5}'.format('vmax'), end=' ')
                    print('{0:^5}'.format('vmin'), end=' ')
                    print('{0:^8}'.format('gvdev'))
                else:
                    print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                    print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                    print('{0:^8.1e}'.format(np.average(net.gen_v_dev)))
            return info_printer

        # NR-based
        ##########
        elif optsolver_name == 'nr':

            def info_printer(solver,header):
                net = solver.problem.wrapped_problem.network
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

            if bus.is_regulated_by_shunt() and not bus.is_slack():
                
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

            if bus.is_regulated_by_tran() and not bus.is_slack():
                
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
