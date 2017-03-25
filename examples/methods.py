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

method = gridopt.power_flow.new_method('NRPF')

method.set_parameters({'quiet': True, 'feastol': 1e-4})

method.solve(net)

results = method.get_results()

print((results['status']))

print((results['iterations']))

problem = results['problem']
problem.show()

print((results['net properties']['bus_v_max']))

method.update_network(net)

print(('%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)))
