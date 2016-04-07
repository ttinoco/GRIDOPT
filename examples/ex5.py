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

gen = net.get_gen(0)
branch = net.get_branch(0)

c1 = pfnet.Contingency(gens=[gen])
c2 = pfnet.Contingency(branches=[branch])

method = gridopt.power_flow.new_method('DCOPF_Corr')
method.solve(net,[c1,c2])
method.update_network(net)

print net.gen_P_cost

