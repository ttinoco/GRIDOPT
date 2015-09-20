#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from nr import NRPF
from dc import DCPF
from augl import AugLPF

def new_method(name):
    
    if name == 'NRPF':
        return NRPF()
    elif name == 'DCPF':
        return DCPF()
    elif name == 'AugLPF':
        return AugLPF()
    else:
        raise ValueError('invalid PF method name')
        
    
    
