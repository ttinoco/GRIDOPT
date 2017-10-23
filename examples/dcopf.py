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

method = gridopt.power_flow.new_method('DCOPF')

method.solve(net)

print((method.results['solver status']))

method.update_network(net)

# generation cost ($/hour)
print((net.gen_P_cost))

# sensitivity to power balance 
bus = net.get_bus(4)
print(("bus %2d %.2e" %(bus.index,bus.sens_P_balance)))

# sensitivity to flow limits 
branch = net.get_branch(6)
print(("branch %2d %.2e %.2e" %(branch.index,
                               branch.sens_P_u_bound,
                               branch.sens_P_l_bound)))

# sensitivity to  gen limits
gen = net.get_generator(2)
print(("gen %2d %.2e %.2e" %(gen.index,
                            gen.sens_P_u_bound,
                            gen.sens_P_l_bound)))
