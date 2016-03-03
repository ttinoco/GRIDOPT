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

    def tearDown(self):
        
        pass

