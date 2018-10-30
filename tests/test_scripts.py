#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import os
import pfnet
import unittest
import numpy as np
import gridopt as gopt

class TestScripts(unittest.TestCase):
    
    def setUp(self):
                    
        pass

    def test_gridopt(self):

        try:
            
            case = os.path.join('.', 'tests', 'resources', 'cases', 'ieee300.mat')
            
            net = pfnet.Parser(case).parse(case)
            
            method = gopt.power_flow.ACPF()
            method.set_parameters({'quiet': True, 'solver': 'nr'})
            method.solve(net)
            method.update_network(net)
            
            gopt.scripts.gridopt.main([case, 'ACPF', '--parameters', 'solver=nr', '--write', 'foo.json', '--quiet'])

            p = pfnet.ParserJSON()
            net_json = p.parse('foo.json')
            pfnet.tests.utils.compare_networks(self, net, net_json)

            gopt.scripts.gridopt.main([case, 'ACPF', '--parameters', 'solver=nr', '--write', 'foo.m', '--quiet'])

            p = pfnet.PyParserMAT()
            net_m = p.parse('foo.m')
            pfnet.tests.utils.compare_networks(self, net, net_m)
            
        finally:
            if os.path.isfile('foo.json'):
                os.remove('foo.json')
            if os.path.isfile('foo.m'):
                os.remove('foo.m')
            
        
