#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import pfnet
import numpy as np
from .method_error import *
from .method import PFmethod
from numpy.linalg import norm
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
        parameters = OptSolverNR.parameters.copy() # solver parameters
        parameters.update(NRPF.parameters) # method parameters
        self.parameters = parameters
        
    def apply_shunt_v_regulation(self,solver):

        # Local variables
        dsus = self.parameters['dsus']
        step = self.parameters['shunt_step']
        p = solver.problem
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
                    vmax = bus.v_max
                    vmin = bus.v_min

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
                    vmax = bus.v_max
                    vmin = bus.v_min
                    
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
                                    net.get_num_tap_changers() +
                                    net.get_num_switched_shunts())*net.num_periods)
            assert(net.num_fixed == (net.get_num_buses_reg_by_gen() +
                                     net.get_num_tap_changers() +
                                     net.get_num_switched_shunts())*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)

        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.Constraint('AC power balance',net))
        problem.add_constraint(pfnet.Constraint('generator active power participation',net))
        problem.add_constraint(pfnet.Constraint('generator reactive power participation',net))
        problem.add_constraint(pfnet.Constraint('variable fixing',net))
        if limit_gens:
            problem.add_heuristic(pfnet.HEUR_TYPE_PVPQ)
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
                norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_tran_v_regulation(s)
            
        def c2(s):
            if (s.k != 0 and
                (not lock_shunts) and 
                norm(s.problem.f,np.inf) < 100.*feastol):
                self.apply_shunt_v_regulation(s)

        def c3(s):
            if s.k > 0:
                s.problem.apply_heuristics(s.x)

        # Termination
        def t1(s):
            if np.min(s.problem.network.bus_v_min) < vmin_thresh:
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

            # Restore net properties
            net.update_properties()
    
    def update_network(self,net):
        
        # Get data
        problem = self.results['problem']
        x = self.results['primal variables']
        lam,nu,mu,pi = self.results['dual variables']
       
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == x.shape)
        assert(net.num_vars == x.size)
        assert(lam is None or not lam.size)
        assert(nu is None or not nu.size)
        assert(mu is None or not mu.size)
        assert(pi is None or not pi.size)

        # Network quantities
        net.set_var_values(x)

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()

