#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import gridopt

net = pfnet.ParserMAT().parse('../tests/resources/cases/ieee14.mat')

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

method = gridopt.power_flow.new_method('ACPF')

method.set_parameters({'quiet': True})

method.solve(net)

results = method.get_results()

print((results['status']))

method.update_network(net)

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))

