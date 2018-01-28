#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import copy
import unittest
import numpy as np
import pfnet as pf
from . import utils
import gridopt as gopt
from numpy.linalg import norm

class TestPowerFlow(unittest.TestCase):
    
    def setUp(self):
                    
        pass

    def test_ACOPF_parameters(self):

        acopf = gopt.power_flow.ACOPF()

        # Set parameters
        acopf.set_parameters({'weight_cost': 1.1,
                              'weight_vmag' : 1.2,
                              'weight_vang' : 1.3,
                              'weight_pq' : 1.4,
                              'weight_t' : 1.5,
                              'weight_b' : 1.6,
                              'thermal_limits' : True,
                              'vmin_thresh' : 0.123,
                              'solver' : 'inlp',
                              'beta_small' : 8.9,
                              'hessian_approximation' : 'test',
                              'maxiter': 432})

        # Check exception for invalid param
        self.assertRaises(gopt.power_flow.method_error.PFmethodError_BadParams,
                          acopf.set_parameters,{'foo' : 'bar'})

        # Get parameters
        params = acopf.get_parameters()

        # Check that set_parameters worked
        self.assertEqual(params['weight_cost'],1.1)
        self.assertEqual(params['weight_vmag'],1.2)
        self.assertEqual(params['weight_vang'],1.3)
        self.assertEqual(params['weight_pq'],1.4)
        self.assertEqual(params['weight_t'],1.5)
        self.assertEqual(params['weight_b'],1.6)
        self.assertEqual(params['thermal_limits'],True)
        self.assertEqual(params['vmin_thresh'],0.123)
        self.assertEqual(params['solver'],'inlp')
        self.assertEqual(params['solver_parameters']['augl']['beta_small'],8.9)
        self.assertEqual(params['solver_parameters']['augl']['maxiter'],432)
        self.assertEqual(params['solver_parameters']['inlp']['maxiter'],432)
        self.assertEqual(params['solver_parameters']['ipopt']['hessian_approximation'],'test')

        # Make a deep copy
        new_params = copy.deepcopy(params)

        # Set manually solver parametres
        new_params['solver_parameters']['augl']['beta_small'] = 1.333e-5
        new_params['solver_parameters']['inlp']['maxiter'] = 555
        new_params['solver_parameters']['ipopt']['hessian_approximation'] = 'new test'

        # Check that this is a separate params dictionary
        self.assertNotEqual(new_params['solver_parameters']['augl']['beta_small'],
                            params['solver_parameters']['augl']['beta_small'])
        self.assertNotEqual(new_params['solver_parameters']['inlp']['maxiter'],
                            params['solver_parameters']['inlp']['maxiter'])
        self.assertNotEqual(new_params['solver_parameters']['ipopt']['hessian_approximation'],
                            params['solver_parameters']['ipopt']['hessian_approximation'])

        # Test setting parameters that specify solver parameters
        acopf.set_parameters(new_params)

        # Check that setting solver parameters worked
        self.assertEqual(new_params['solver_parameters']['augl']['beta_small'],
                         params['solver_parameters']['augl']['beta_small'])
        self.assertEqual(new_params['solver_parameters']['inlp']['maxiter'],
                         params['solver_parameters']['inlp']['maxiter'])
        self.assertEqual(new_params['solver_parameters']['ipopt']['hessian_approximation'],
                         params['solver_parameters']['ipopt']['hessian_approximation'])

    def test_DCPF(self):

        for case in utils.test_cases:

            method = gopt.power_flow.new_method('DCPF')
            self.assertTrue(isinstance(method, gopt.power_flow.DCPF))

            method.set_parameters({'solver' : 'superlu'})
            self.assertTrue('solver_parameters' in method.get_parameters().keys())

            net = pf.Parser(case).parse(case)

            method.solve(net)

            results = method.get_results()

            self.assertEqual(results['solver status'], 'solved')
            self.assertTrue(results['solver name'] in ['mumps','superlu'])
            self.assertTrue(isinstance(results['network snapshot'], pf.Network))

    def test_ACPF_nr_heuristics(self):

        for case in utils.test_cases:

            net = pf.Parser(case).parse(case)

            for Q_par in ['range', 'fraction']:

                method = gopt.power_flow.new_method('ACPF')
                method.set_parameters(params={'solver': 'nr',
                                              'quiet': True})
                method.solve(net)

                results = method.get_results()

                net_snap = results['network snapshot']
                self.assertLess(np.abs(net_snap.bus_P_mis), 1e-2) # MW
                self.assertLess(np.abs(net_snap.bus_Q_mis), 1e-2) # MVAr

                eps = 1e-4
                for bus in net_snap.buses:
                    if bus.is_regulated_by_gen() and not bus.is_slack():
                        for gen in bus.reg_generators:
                            self.assertLessEqual(gen.Q, gen.Q_max+1e-10)
                            self.assertGreaterEqual(gen.Q, gen.Q_min-1e-10)
                        if np.abs(bus.v_mag-bus.v_set) < eps: # v at set
                            Qtotal = 0.
                            norm = 0.
                            for gen in bus.reg_generators:
                                if np.abs(gen.Q-gen.Q_max) > eps and np.abs(gen.Q-gen.Q_min) > eps:
                                    Qtotal += gen.Q
                                    norm += gen.Q_par
                            for gen in bus.reg_generators:
                                if np.abs(gen.Q-gen.Q_max) > eps and np.abs(gen.Q-gen.Q_min) > eps:
                                    self.assertLess(np.abs(gen.Q-gen.Q_par*Qtotal/norm), eps)
                        else: # v not at set
                            num = 0
                            for gen in bus.reg_generators:
                                if np.abs(gen.Q-gen.Q_max) <= eps or np.abs(gen.Q-gen.Q_min) <= eps:
                                    num += 1
                            self.assertGreaterEqual(num, 1)                                
                                
    def test_ACPF_solutions(self):

        print('')
        
        T = 2

        sol_types = {'sol1': 'no--controls',
                     'sol2': 'gen-controls',
                     'sol3': 'all-controls'}
        
        for case in utils.test_cases:
            for sol in list(sol_types.keys()):
                for solver in ['nr','augl','ipopt', 'inlp']:
                    
                    method = gopt.power_flow.new_method('ACPF')
                    method.set_parameters(params={'solver': solver})
                    
                    net = pf.Parser(case).parse(case)
                    netMP = pf.Parser(case).parse(case,T)
                    
                    self.assertEqual(net.num_periods,1)
                    self.assertEqual(netMP.num_periods,T)

                    # Only small
                    if net.num_buses > 3000:
                        continue
                    
                    sol_file = utils.get_pf_solution_file(case,sol)
                    sol_data = utils.read_pf_solution_file(sol_file)
                                          
                    # Set parameters
                    if sol == 'sol1':
                        method.set_parameters({'limit_gens': False})
                    elif sol == 'sol2':
                        pass # defaults
                    elif sol == 'sol3':
                        method.set_parameters({'lock_taps': False})
                    else:
                        raise ValueError('invalid solution type')
                    method.set_parameters({'quiet': True})

                    bus_P_mis = net.bus_P_mis
                    try:
                        method.solve(net)
                    except ImportError:
                        continue # no ipopt
                    results = method.get_results()
                    self.assertEqual(results['solver name'], solver)
                    self.assertEqual(results['solver status'],'solved')
                    self.assertEqual(net.bus_P_mis,bus_P_mis)
                    self.assertLessEqual(results['network snapshot'].bus_P_mis,bus_P_mis)
                    method.update_network(net)

                    method.solve(netMP)
                    resultsMP = method.get_results()
                    self.assertEqual(resultsMP['solver status'],'solved')
                    method.update_network(netMP)

                    self.assertLess(norm(resultsMP['network snapshot'].bus_P_mis-netMP.bus_P_mis,np.inf),1e-10)
                    self.assertLess(norm(resultsMP['network snapshot'].bus_Q_mis-netMP.bus_Q_mis,np.inf),1e-10)
                    self.assertLess(norm(resultsMP['network snapshot'].gen_P_cost-netMP.gen_P_cost,np.inf),1e-10)

                    print("%s\t%s\t%s\t%d" %(case.split('/')[-1],
                                             sol_types[sol],
                                             solver,
                                             results['solver iterations']))

                    # No sol
                    if sol_data is None:
                        continue

                    v_mag_tol = sol_data['v_mag_tol']
                    v_ang_tol = sol_data['v_ang_tol']
                    bus_data = sol_data['bus_data']

                    counter = 0
                    v_mag_error = []
                    v_ang_error = []
                    for bus_num,val in list(bus_data.items()):
                        
                        v_mag = val['v_mag']
                        v_ang = val['v_ang']
                        
                        try:
                            busMP = netMP.get_bus_from_number(bus_num)
                            bus = net.get_bus_from_number(bus_num)
                        except pf.NetworkError:
                            continue
                                                       
                        for t in range(T):
                            v_mag_error.append(np.abs(busMP.v_mag[t]-v_mag))
                            v_ang_error.append(np.abs(busMP.v_ang[t]*180./np.pi-v_ang))
                        v_mag_error.append(np.abs(bus.v_mag-v_mag))
                        v_ang_error.append(np.abs(bus.v_ang*180./np.pi-v_ang))
                       
                        counter += 1

                    self.assertEqual(len(v_mag_error),len(v_ang_error))
                    if len(v_mag_error) > 0:
                        self.assertLessEqual(np.max(v_mag_error),v_mag_tol)
                        self.assertLessEqual(np.max(v_ang_error),v_ang_tol)
                    self.assertEqual(len(v_mag_error),counter*(T+1))
                    self.assertEqual(len(v_ang_error),counter*(T+1))

    def test_ACOPF_solutions(self):

        print('')
        
        eps = 0.5 # %
        
        method_ipopt = gopt.power_flow.new_method('ACOPF')
        method_ipopt.set_parameters(params={'solver':'ipopt','quiet': True})
        method_augl = gopt.power_flow.new_method('ACOPF')
        method_augl.set_parameters(params={'solver':'augl','quiet': True})
        method_inlp = gopt.power_flow.new_method('ACOPF')
        method_inlp.set_parameters(params={'solver':'inlp','quiet': True})

        skipcases = ['aesoSL2014.raw','case3012wp.mat','case9241.mat','case32.art']
            
        for case in utils.test_cases:

            if case.split('/')[-1] in skipcases:
                continue

            net = pf.Parser(case).parse(case)

            # Only small
            if net.num_buses > 3000:
                continue
            
            self.assertEqual(net.num_periods,1)

            # IPOPT
            try:
                net.update_properties()
                gen_P_cost = net.gen_P_cost
                method_ipopt.solve(net)
                has_ipopt = True
                self.assertEqual(method_ipopt.results['solver status'],'solved')
                self.assertEqual(net.gen_P_cost,gen_P_cost)
                self.assertNotEqual(method_ipopt.results['network snapshot'].gen_P_cost,gen_P_cost)
                self.assertEqual(method_ipopt.results['solver name'], 'ipopt')
                x1 = method_ipopt.get_results()['solver primal variables']
                i1 = method_ipopt.get_results()['solver iterations']
                p1 = method_ipopt.get_results()['network snapshot'].gen_P_cost
                print("%s\t%s\t%d" %(case.split('/')[-1],'ipopt',i1))
            except ImportError:
                has_ipopt = False

            # INLP
            net.update_properties()
            gen_P_cost = net.gen_P_cost
            method_inlp.solve(net)
            self.assertEqual(method_inlp.results['solver status'],'solved')
            self.assertEqual(net.gen_P_cost,gen_P_cost)
            self.assertNotEqual(method_inlp.results['network snapshot'].gen_P_cost,gen_P_cost)
            self.assertEqual(method_inlp.results['solver name'], 'inlp')
            x2 = method_inlp.get_results()['solver primal variables']
            i2 = method_inlp.get_results()['solver iterations']
            p2 = method_inlp.get_results()['network snapshot'].gen_P_cost
            print("%s\t%s\t%d" %(case.split('/')[-1],'inlp',i2))

            # AUGL
            net.update_properties()
            gen_P_cost = net.gen_P_cost
            method_augl.solve(net)
            self.assertEqual(method_augl.results['solver status'],'solved')
            self.assertEqual(net.gen_P_cost,gen_P_cost)
            self.assertNotEqual(method_augl.results['network snapshot'].gen_P_cost,gen_P_cost)
            self.assertEqual(method_augl.results['solver name'], 'augl')
            x3 = method_augl.get_results()['solver primal variables']
            i3 = method_augl.get_results()['solver iterations']
            p3 = method_augl.get_results()['network snapshot'].gen_P_cost
            print("%s\t%s\t%d" %(case.split('/')[-1],'augl',i3))

            # Checks
            if has_ipopt:
                error = 100*(p1-p3)/abs(p3)
                self.assertLess(np.abs(error),eps)
                self.assertNotEqual(p1,p3)
            error = 100*(p2-p3)/abs(p3)
            self.assertLess(np.abs(error),eps)
            self.assertNotEqual(p2,p3)

    def test_DCOPF_solutions(self):

        T = 2

        infcases = ['ieee25.raw']

        skipcases = ['case1354.mat','case2869.mat',
                     'case3375wp.mat','case9241.mat']
        
        method = gopt.power_flow.new_method('DCOPF')

        for case in utils.test_cases:

            if case.split('/')[-1] in skipcases:
                continue
        
            net = pf.Parser(case).parse(case,T)

            for branch in net.branches:
                if branch.ratingA == 0:
                    branch.ratingA = 100
            
            self.assertEqual(net.num_periods,T)
            
            method.set_parameters({'quiet':True, 
                                   'tol': 1e-6,
                                   'thermal_limits': True})

            try:
                net.update_properties()
                gen_P_cost = net.gen_P_cost
                method.solve(net)
                self.assertEqual(method.results['solver status'],'solved')
                self.assertEqual(method.results['solver name'], 'iqp')
                self.assertTrue(np.all(net.gen_P_cost == gen_P_cost))
                self.assertTrue(np.all(method.results['network snapshot'].gen_P_cost != gen_P_cost))
            except gopt.power_flow.PFmethodError:
                self.assertTrue(case.split('/')[-1] in infcases)
                self.assertEqual(method.results['solver status'],'error')
                
            results = method.get_results()
                
            method.update_network(net)
           
            self.assertLess(norm(results['network snapshot'].bus_P_mis-net.bus_P_mis,np.inf),1e-10)
            self.assertLess(norm(results['network snapshot'].bus_Q_mis-net.bus_Q_mis,np.inf),1e-10)
            self.assertLess(norm(results['network snapshot'].gen_P_cost-net.gen_P_cost,np.inf),1e-10)

            gen_P_cost0 = net.gen_P_cost
            load_P_util0 = net.load_P_util
            self.assertTupleEqual(gen_P_cost0.shape,(T,))
            self.assertTupleEqual(load_P_util0.shape,(T,))
            
            x = results['solver primal variables']
            lam0,nu0,mu0,pi0 = results['solver dual variables']

            self.assertTupleEqual(x.shape,((net.num_buses-
                                            net.get_num_slack_buses()+
                                            net.get_num_P_adjust_gens())*T,))
            self.assertTupleEqual(x.shape,(net.num_vars,))
            self.assertTupleEqual(lam0.shape,(net.num_buses*T,))
            self.assertTrue(nu0.size == 0)
            self.assertTupleEqual(mu0.shape,(net.num_vars+net.num_branches*T,))
            self.assertTupleEqual(pi0.shape,(net.num_vars+net.num_branches*T,))

            # Network update (vars and sensitivities)
            xx = x[:net.num_vars]
            for t in range(T):
                for bus in net.buses:
                    if not bus.is_slack():
                        self.assertEqual(bus.v_ang[t],xx[bus.index_v_ang[t]])
                        self.assertEqual(bus.sens_v_ang_u_bound[t],mu0[bus.index_v_ang[t]])
                        self.assertEqual(bus.sens_v_ang_l_bound[t],pi0[bus.index_v_ang[t]])
                for gen in net.generators:
                    if gen.is_P_adjustable():
                        self.assertEqual(gen.P[t],xx[gen.index_P[t]])
                        self.assertEqual(gen.sens_P_u_bound[t],mu0[gen.index_P[t]])
                        self.assertEqual(gen.sens_P_l_bound[t],pi0[gen.index_P[t]])
                for branch in net.branches:
                    self.assertEqual(branch.sens_P_u_bound[t],mu0[net.num_vars+branch.index+t*net.num_branches])
                    self.assertEqual(branch.sens_P_l_bound[t],pi0[net.num_vars+branch.index+t*net.num_branches])

            # No thermal limits
            method.set_parameters({'thermal_limits':False})
            method.solve(net)
            self.assertEqual(method.results['solver status'],'solved')
            results = method.get_results()
            method.update_network(net)
            gen_P_cost1 = net.gen_P_cost
            load_P_util1 = net.load_P_util
            lam1,nu1,mu1,pi1 = results['solver dual variables']
            self.assertTupleEqual(mu1.shape,(net.num_vars,))
            self.assertTupleEqual(pi1.shape,(net.num_vars,))

            # Elastic loads
            for load in net.loads:
                load.P_max = load.P+1.
                load.P_min = load.P-1.
            self.assertEqual(net.get_num_P_adjust_loads(),net.num_loads)
            for load in net.loads:
                self.assertFalse(load.has_flags('variable','active power'))
                self.assertFalse(load.has_flags('bounded','active power'))
                self.assertTrue(np.all(load.P_min < load.P_max))
            method.solve(net)
            for load in net.loads:
                self.assertFalse(load.has_flags('variable','active power'))
                self.assertFalse(load.has_flags('bounded','active power'))
            results = method.get_results()
            method.update_network(net)
            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertTrue(load.has_flags('bounded','active power'))
            self.assertEqual(method.results['solver status'],'solved')
            self.assertTrue(np.all(net.gen_P_cost-net.load_P_util < gen_P_cost1-load_P_util1))

            x = results['solver primal variables']
            lam2,nu2,mu2,pi2 = results['solver dual variables']

            self.assertTupleEqual(x.shape,((net.get_num_P_adjust_loads()+
                                            net.num_buses-
                                            net.get_num_slack_buses()+
                                            net.get_num_P_adjust_gens())*net.num_periods,))
            self.assertTupleEqual(x.shape,(net.num_vars,))
            self.assertTupleEqual(lam2.shape,(net.num_buses*net.num_periods,))
            self.assertTrue(nu2.size == 0)
            self.assertTupleEqual(mu2.shape,x.shape)
            self.assertTupleEqual(pi2.shape,x.shape)

            xx = x[:net.num_vars]
            for t in range(T):
                for bus in net.buses:
                    if not bus.is_slack():
                        self.assertEqual(bus.v_ang[t],xx[bus.index_v_ang[t]])
                        self.assertEqual(bus.sens_v_ang_u_bound[t],mu2[bus.index_v_ang[t]])
                        self.assertEqual(bus.sens_v_ang_l_bound[t],pi2[bus.index_v_ang[t]])
                for gen in net.generators:
                    if gen.is_P_adjustable():
                        self.assertEqual(gen.P[t],xx[gen.index_P[t]])
                        self.assertEqual(gen.sens_P_u_bound[t],mu2[gen.index_P[t]])
                        self.assertEqual(gen.sens_P_l_bound[t],pi2[gen.index_P[t]])
                for load in net.loads:
                    if load.is_P_adjustable():
                        self.assertEqual(load.P[t],xx[load.index_P[t]])
                        self.assertEqual(load.sens_P_u_bound[t],mu2[load.index_P[t]])
                        self.assertEqual(load.sens_P_l_bound[t],pi2[load.index_P[t]])
                     
    def tearDown(self):
        
        pass

