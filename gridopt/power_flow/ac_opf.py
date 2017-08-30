#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import numpy as np
from .method_error import *
from .method import PFmethod
        
class ACOPF(PFmethod):
    """
    AC optimal power flow method.
    """

    name = 'ACOPF'

    parameters = {'weight_cost' : 1e0,     # weight for generation cost
                  'weight_vmag' : 0.,      # weight for voltage magnitude regularization
                  'weight_vang' : 0.,      # weight for voltage angle regularization
                  'weight_pq' : 0.,        # weight for generator power regularization
                  'weight_t' : 0.,         # weight for tap ratios regularization
                  'weight_b' : 0.,         # weight for shunt susceptances regularization
                  'thermal_limits': False, # flag for thermal limits
                  'vmin_thresh': 0.1,      # threshold for vmin termination
                  'optsolver': 'augl'}     # OPTALG optimization solver (augl,ipopt)

    parameters_augl = {'feastol' : 1e-4,
                       'optol' : 1e-4,
                       'kappa' : 1e-2}

    parameters_ipopt = {}

    parameters_inlp = {}

    def __init__(self):

        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverINLP

        # Parent init
        PFmethod.__init__(self)

        # Optsolver params
        augl_params = OptSolverAugL.parameters.copy()
        augl_params.update(self.parameters_augl)   # overwrite defaults
        ipopt_params = OptSolverIpopt.parameters.copy()
        ipopt_params.update(self.parameters_ipopt) # overwrite defaults
        inlp_params = OptSolverINLP.parameters.copy()
        inlp_params.update(self.parameters_inlp) # overwrite defaults
        self.parameters.update(ACOPF.parameters)
        self.parameters['optsolver_params'] = {'augl': augl_params,
                                               'ipopt': ipopt_params,
                                               'inlp': inlp_params}
                   
    def create_problem(self,net):
        
        import pfnet

        # Parameters
        params = self.parameters
        wcost  = params['weight_cost']
        wvmag  = params['weight_vmag']
        wvang = params['weight_vang']
        wpq = params['weight_pq']
        wt = params['weight_t']
        wb = params['weight_b']        
        th = params['thermal_limits']
        
        # Clear flags
        net.clear_flags()
        
        # Voltage magnitudes
        net.set_flags('bus',
                      ['variable','bounded'], 
                      'any',
                      'voltage magnitude')
        
        # Voltage angles
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')

        # Generator powers
        net.set_flags('generator',
                      ['variable','bounded'],
                      'not on outage',
                      ['active power','reactive power'])

        try:
            assert(net.num_vars == (2*net.num_buses-net.get_num_slack_buses() +
                                    2*net.get_num_gens_not_on_outage())*net.num_periods)
            assert(net.num_bounded == (2*net.get_num_gens_not_on_outage() +
                                       net.num_buses)*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
                                    
        # Problem
        problem = pfnet.Problem(net)

        # Constraints
        problem.add_constraint(pfnet.Constraint('AC power balance',net))
        problem.add_constraint(pfnet.Constraint('variable bounds',net))
        if th:
            problem.add_constraint(pfnet.Constraint("AC branch flow limits",net))

        # Functions
        problem.add_function(pfnet.Function('generation cost',
                                            wcost/max([net.num_generators,1.]),net))
        if wvmag:
            problem.add_function(pfnet.Function('soft voltage magnitude limits',
                                                wvmag/max([net.num_buses,1.]),net))
        if wvang:
            problem.add_function(pfnet.Function('voltage angle regularization',
                                                wvang/max([net.num_buses,1.]),net))
        if wpq:
            problem.add_function(pfnet.Function('generator powers regularization',
                                                wpq/max([net.num_generators,1.]),net))
        if wt:
            problem.add_function(pfnet.Function('tap ratio regularization',
                                                wt/max([net.get_num_tap_changers(),1.]),net))
        if wb:
            problem.add_function(pfnet.Function('susceptance regularization',
                                                wb/max([net.get_num_switched_shunts(),1.]),net))
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net):

        from optalg.opt_solver import OptSolverError, OptTermination
        from optalg.opt_solver import OptSolverAugL, OptSolverIpopt, OptSolverINLP
        
        # Parameters
        params = self.parameters
        vmin_thresh = params['vmin_thresh']
        optsolver_name = params['optsolver']
        optsolver_params = params['optsolver_params']

        # Opt solver
        if optsolver_name == 'augl':
            optsolver = OptSolverAugL()
        elif optsolver_name == 'ipopt':
            optsolver = OptSolverIpopt()
        elif optsolver_name == 'inlp':
            optsolver = OptSolverINLP()
        else:
            raise PFmethodError_BadOptSolver(self)
        optsolver.set_parameters(optsolver_params[optsolver_name])
        
        # Problem
        problem = self.create_problem(net)

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
        try:
            optsolver.solve(problem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:
            
            # Get results
            self.set_status(optsolver.get_status())
            self.set_error_msg(optsolver.get_error_msg())
            self.set_iterations(optsolver.get_iterations())
            self.set_primal_variables(optsolver.get_primal_variables())
            self.set_dual_variables(optsolver.get_dual_variables())
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)

            # Restors net properties
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
        assert(net.num_vars+problem.num_extra_vars == x.size)
        assert(problem.A.shape[0] == lam.size)
        assert(problem.f.shape[0] == nu.size)
        assert(problem.G.shape[0] == mu.size)
        assert(problem.G.shape[0] == pi.size)

        # Network quantities
        net.set_var_values(x[:net.num_vars])

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam,nu,mu,pi)

    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.wrapped_problem.network
            if header:
                print('{0:^5}'.format('vmax'), end=' ')
                print('{0:^5}'.format('vmin'), end=' ')
                print('{0:^6}'.format('bvvio'), end=' ')
                print('{0:^6}'.format('gQvio'), end=' ')
                print('{0:^6}'.format('gPvio'))
            else:
                print('{0:^5.2f}'.format(np.average(net.bus_v_max)), end=' ')
                print('{0:^5.2f}'.format(np.average(net.bus_v_min)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.bus_v_vio)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.gen_Q_vio)), end=' ')
                print('{0:^6.0e}'.format(np.average(net.gen_P_vio)))
        return info_printer

