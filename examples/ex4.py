#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import sys
import pfnet
import numpy as np

net = pfnet.Network()
net.load(sys.argv[1])

T = 3

load_data = {}
for load in net.loads:
    load_data[load.index] = np.random.rand(T)

def net_modifier(net,t):
    print 'modifying net for time %d' %t
    for load in net.loads:
        load.P = load_data[load.index][t]
        load.P_max = 1.05*load.P
        load.P_min = 0.95*load.P

map(lambda t: net_modifier(net,t),range(T))


