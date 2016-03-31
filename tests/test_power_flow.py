#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import utils
import unittest
import numpy as np
import pfnet as pf
import gridopt as gopt

INFCASE = './tests/resources/ieee25.raw'

class TestPowerFlow(unittest.TestCase):
    
    def setUp(self):
        
        # Network
        self.net = pf.Network()

    def test_method_solutions(self):

        print ''
        net = self.net
        sol_types = {'sol1': 'no controls',
                     'sol2': 'gen controls',
                     'sol3': 'all controls'}
        
        for method_name in ['NRPF','AugLPF']:
            for case in utils.test_cases:
                for sol in sol_types.keys():

                    method = gopt.power_flow.new_method(method_name)
                    
                    net.load(case)
                    
                    sol_data = utils.read_solution_data(case+'.'+sol)
                
                    # Skip
                    if (sol == 'sol1' and         # no controls
                        method_name == 'AugLPF'):

                        continue

                    if sol == 'sol1':   # no controls
                        method.set_parameters({'limit_gens': False})
                    elif sol == 'sol2': # generator voltage control
                        pass
                    else:               # generator and tap controls
                        method.set_parameters({'lock_taps': False})
                    method.set_parameters({'quiet': True})

                    method.solve(net)

                    results = method.get_results()

                    self.assertEqual(results['status'],'solved')
                    
                    method.update_network(net)

                    self.assertLess(np.abs(results['net_properties']['bus_P_mis']-net.bus_P_mis),1e-10)
                    self.assertLess(np.abs(results['net_properties']['bus_Q_mis']-net.bus_Q_mis),1e-10)
                    self.assertLess(np.abs(results['net_properties']['gen_P_cost']-net.gen_P_cost),1e-10)

                    v_mag_tol = sol_data['v_mag_tol']
                    v_ang_tol = sol_data['v_ang_tol']
                    bus_data = sol_data['bus_data']

                    v_mag_error = [0]
                    v_ang_error = [0]
                    for bus_num,val in bus_data.items():
                        
                        v_mag = val['v_mag']
                        v_ang = val['v_ang']
                        
                        try:
                            bus = net.get_bus_by_number(bus_num)
                        except pf.NetworkError:
                            continue

                        v_mag_error.append(np.abs(bus.v_mag-v_mag))
                        v_ang_error.append(np.abs(bus.v_ang*180./np.pi-v_ang))
                    
                    print method_name,case,sol_types[sol],len(v_mag_error),len(v_ang_error)

                    self.assertLessEqual(np.max(v_mag_error),v_mag_tol)
                    self.assertLessEqual(np.max(v_ang_error),v_ang_tol)

    def test_DCOPF(self):
        
        net = self.net
        method = gopt.power_flow.new_method('DCOPF')

        for case in utils.test_cases:
        
            net.load(case)

            method.set_parameters({'quiet':True})

            try:
                method.solve(net)
                self.assertEqual(method.results['status'],'solved')
            except gopt.power_flow.PFmethodError:
                self.assertEqual(case,INFCASE)
                self.assertEqual(method.results['status'],'error')
                
            results = method.get_results()
                
            method.update_network(net)
           
            self.assertLess(np.abs(results['net_properties']['bus_P_mis']-net.bus_P_mis),1e-10)
            self.assertLess(np.abs(results['net_properties']['bus_Q_mis']-net.bus_Q_mis),1e-10)
            self.assertLess(np.abs(results['net_properties']['gen_P_cost']-net.gen_P_cost),1e-10)
            
            x = results['primal_variables']
            lam,nu,mu,pi = results['dual_variables']

            self.assertTupleEqual(x.shape,(net.num_branches+
                                           net.num_buses-
                                           net.get_num_slack_buses()+
                                           net.get_num_P_adjust_gens(),))
            self.assertTupleEqual(x.shape,(net.num_vars+net.num_branches,))
            self.assertTupleEqual(lam.shape,(net.num_buses+net.num_branches,))
            self.assertTrue(nu is None)
            self.assertTupleEqual(mu.shape,x.shape)
            self.assertTupleEqual(pi.shape,x.shape)            

            xx = x[:net.num_vars]
            for bus in net.buses:
                if not bus.is_slack():
                    self.assertEqual(bus.v_ang,xx[bus.index_v_ang])
                    self.assertEqual(bus.sens_v_ang_u_bound,mu[bus.index_v_ang])
                    self.assertEqual(bus.sens_v_ang_l_bound,pi[bus.index_v_ang])
            for gen in net.generators:
                if gen.is_P_adjustable():
                    self.assertEqual(gen.P,xx[gen.index_P])
                    self.assertEqual(gen.sens_P_u_bound,mu[gen.index_P])
                    self.assertEqual(gen.sens_P_l_bound,pi[gen.index_P])

    def test_DCOPF_prev(self):
        
        net = self.net
        method = gopt.power_flow.new_method('PreventiveDCOPF')

        for case in utils.test_cases:
        
            net.load(case)
            
            cont = pf.Contingency()
            cont.add_gen_outage(net.get_gen(0))

            method.set_parameters({'quiet':True})

            try:
                method.solve(net,[cont])
            except Exception:
                raise

    def tearDown(self):
        
        pass

