#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
from . import utils
import unittest
import numpy as np
import pfnet as pf
import gridopt as gopt
from numpy.linalg import norm

class TestPowerFlow(unittest.TestCase):
    
    def setUp(self):
                    
        pass

    def test_ACPF(self):

        print('')
        
        T = 2

        sol_types = {'sol1': 'no--controls',
                     'sol2': 'gen-controls',
                     'sol3': 'all-controls'}
        
        for case in utils.test_cases:
            for sol in list(sol_types.keys()):
                for optsolver in ['nr','augl','ipopt']:
                    
                    method = gopt.power_flow.new_method('ACPF')
                    method.set_parameters(params={'optsolver':optsolver})
                    
                    net = pf.Parser(case).parse(case)
                    netMP = pf.Parser(case).parse(case,T)
                    
                    self.assertEqual(net.num_periods,1)
                    self.assertEqual(netMP.num_periods,T)

                    sol_file = utils.get_pf_solution_file(case,sol)
                    sol_data = utils.read_pf_solution_file(sol_file)

                    # No sol
                    if sol_data is None:
                        continue
                                          
                    # Skip (no controls is not supported by augl,ipopt formulation
                    if sol == 'sol1' and optsolver in ['augl','ipopt']:
                        continue

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
                    self.assertEqual(results['status'],'solved')
                    self.assertEqual(net.bus_P_mis,bus_P_mis)
                    self.assertNotEqual(results['net properties']['bus_P_mis'],bus_P_mis)
                    method.update_network(net)

                    method.solve(netMP)
                    resultsMP = method.get_results()
                    self.assertEqual(resultsMP['status'],'solved')
                    method.update_network(netMP)

                    self.assertLess(norm(resultsMP['net properties']['bus_P_mis']-netMP.bus_P_mis,np.inf),1e-10)
                    self.assertLess(norm(resultsMP['net properties']['bus_Q_mis']-netMP.bus_Q_mis,np.inf),1e-10)
                    self.assertLess(norm(resultsMP['net properties']['gen_P_cost']-netMP.gen_P_cost,np.inf),1e-10)

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
                            busMP = netMP.get_bus_by_number(bus_num)
                            bus = net.get_bus_by_number(bus_num)
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
                        print("%s\t%s\t%s" %(optsolver,case.split('/')[-1],sol_types[sol]))
                    self.assertEqual(len(v_mag_error),counter*(T+1))
                    self.assertEqual(len(v_ang_error),counter*(T+1))

    def test_ACOPF(self):

        print('')
        
        eps = 0.5 # %
        
        method_ipopt = gopt.power_flow.new_method('ACOPF')
        method_ipopt.set_parameters(params={'optsolver':'ipopt'})
        method_augl = gopt.power_flow.new_method('ACOPF')
        method_augl.set_parameters(params={'optsolver':'augl'})

        skipcases = ['aesoSL2014.raw','case3012wp.mat','case9241.mat','case32.art']
            
        for case in utils.test_cases:

            if case.split('/')[-1] in skipcases:
                continue

            net = pf.Parser(case).parse(case)
            
            self.assertEqual(net.num_periods,1)
            
            method_ipopt.set_parameters({'quiet':True})
            method_augl.set_parameters({'quiet':True})

            try:
                net.update_properties()
                gen_P_cost = net.gen_P_cost
                method_ipopt.solve(net)
                has_ipopt = True
                self.assertEqual(method_ipopt.results['status'],'solved')
                self.assertEqual(net.gen_P_cost,gen_P_cost)
                self.assertNotEqual(method_ipopt.results['net properties']['gen_P_cost'],gen_P_cost)
                x1 = method_ipopt.get_results()['primal variables']
                i1 = method_ipopt.get_results()['iterations']
                p1 = method_ipopt.get_results()['net properties']['gen_P_cost']
                print("%s\t%s\t%d" %('ipopt',case.split('/')[-1],i1))
            except ImportError:
                has_ipopt = False
            
            net.update_properties()
            gen_P_cost = net.gen_P_cost
            method_augl.solve(net)
            self.assertEqual(method_augl.results['status'],'solved')
            self.assertEqual(net.gen_P_cost,gen_P_cost)
            self.assertNotEqual(method_augl.results['net properties']['gen_P_cost'],gen_P_cost)
            x2 = method_augl.get_results()['primal variables']
            i2 = method_augl.get_results()['iterations']
            p2 = method_augl.get_results()['net properties']['gen_P_cost']
            print("%s\t%s\t%d" %('augl',case.split('/')[-1],i2))
            if has_ipopt:
                error = 100*(p1-p2)/abs(p2)
                self.assertLess(np.abs(error),eps)
                self.assertNotEqual(p1,p2)

    @unittest.skip("temporarily disabled")            
    def test_DCOPF(self):

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
                                   'thermal_factor': 0.93,
                                   'tol': 1e-6})

            try:
                net.update_properties()
                gen_P_cost = net.gen_P_cost
                method.solve(net)
                self.assertEqual(method.results['status'],'solved')
                self.assertTrue(np.all(net.gen_P_cost == gen_P_cost))
                self.assertTrue(np.all(method.results['net properties']['gen_P_cost'] != gen_P_cost))
            except gopt.power_flow.PFmethodError:
                self.assertTrue(case.split('/')[-1] in infcases)
                self.assertEqual(method.results['status'],'error')
                
            results = method.get_results()
                
            method.update_network(net)
           
            self.assertLess(norm(results['net properties']['bus_P_mis']-net.bus_P_mis,np.inf),1e-10)
            self.assertLess(norm(results['net properties']['bus_Q_mis']-net.bus_Q_mis,np.inf),1e-10)
            self.assertLess(norm(results['net properties']['gen_P_cost']-net.gen_P_cost,np.inf),1e-10)

            gen_P_cost0 = net.gen_P_cost
            load_P_util0 = net.load_P_util
            self.assertTupleEqual(gen_P_cost0.shape,(T,))
            self.assertTupleEqual(load_P_util0.shape,(T,))
            
            x = results['primal variables']
            lam0,nu0,mu0,pi0 = results['dual variables']

            self.assertTupleEqual(x.shape,((net.num_branches+
                                            net.num_buses-
                                            net.get_num_slack_buses()+
                                            net.get_num_P_adjust_gens())*T,))
            self.assertTupleEqual(x.shape,(net.num_vars+net.num_branches*T,))
            self.assertTupleEqual(lam0.shape,((net.num_buses+net.num_branches)*T,))
            self.assertTrue(nu0.size == 0)
            self.assertTupleEqual(mu0.shape,x.shape)
            self.assertTupleEqual(pi0.shape,x.shape)

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

            # gen outage 
            if net.get_num_P_adjust_gens() > 1:
                cont = pf.Contingency(gens=[net.get_generator(0)])
                cont.apply()
                try:
                    method.solve(net)
                    self.assertEqual(method.results['status'],'solved')
                except gopt.power_flow.PFmethodError:
                    self.assertTrue(case.split('/')[-1] in infcases)
                    self.assertEqual(method.results['status'],'error')
                cont.clear()

            # no thermal limits
            method.set_parameters({'thermal_limits':False})
            method.solve(net)
            self.assertEqual(method.results['status'],'solved')
            results = method.get_results()
            method.update_network(net)
            gen_P_cost1 = net.gen_P_cost
            load_P_util1 = net.load_P_util
            lam1,nu1,mu1,pi1 = results['dual variables']
            if ((norm(mu0[net.num_vars:],np.inf) > 1e-3 or 
                 norm(pi0[net.num_vars:],np.inf) > 1e-3) and (case.split('/')[-1] not in infcases)):
                self.assertTrue(np.all(gen_P_cost1 <= gen_P_cost0))
            self.assertLess(norm(mu1[net.num_vars:],np.inf),1e-6)
            self.assertLess(norm(pi1[net.num_vars:],np.inf),1e-6)
           
            # elastic loads
            for load in net.loads:
                load.P_max = load.P[0]+1.
                load.P_min = load.P[0]-1.
            for load in net.loads:
                self.assertFalse(load.has_flags('variable','active power'))
                self.assertFalse(load.has_flags('bounded','active power'))
            method.solve(net)
            for load in net.loads:
                self.assertTrue(load.has_flags('variable','active power'))
                self.assertTrue(load.has_flags('bounded','active power'))
                self.assertEqual(method.results['status'],'solved')
            results = method.get_results()
            method.update_network(net)
            self.assertTrue(np.all(net.gen_P_cost-net.load_P_util < gen_P_cost1-load_P_util1))

            x = results['primal variables']
            lam2,nu2,mu2,pi2 = results['dual variables']

            self.assertTupleEqual(x.shape,((net.num_branches+
                                            net.get_num_P_adjust_loads()+
                                            net.num_buses-
                                            net.get_num_slack_buses()+
                                            net.get_num_P_adjust_gens())*net.num_periods,))
            self.assertTupleEqual(x.shape,(net.num_vars+net.num_branches*net.num_periods,))
            self.assertTupleEqual(lam2.shape,((net.num_buses+net.num_branches)*net.num_periods,))
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
                for branch in net.branches:
                    self.assertEqual(branch.sens_P_u_bound[t],mu2[net.num_vars+branch.index+t*net.num_branches])
                    self.assertEqual(branch.sens_P_l_bound[t],pi2[net.num_vars+branch.index+t*net.num_branches])
                     
    def tearDown(self):
        
        pass

