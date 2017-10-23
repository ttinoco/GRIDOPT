#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import sys
sys.path.append('.')

import pfnet
import gridopt

net = pfnet.Parser(sys.argv[1]).parse(sys.argv[1])

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

method = gridopt.power_flow.new_method('ACPF')

method.set_parameters({'quiet': True})

method.solve(net)

results = method.get_results()

print((results['solver status']))

method.update_network(net)

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

