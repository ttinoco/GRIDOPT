#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import sys
import pfnet
import gridopt

net = pfnet.Network()
net.load(sys.argv[1])

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

method = gridopt.power_flow.new_method('NRPF')

method.set_parameters({'quiet': True})

method.solve(net)

results = method.get_results()

print((results['status']))

method.update_network(net)

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

