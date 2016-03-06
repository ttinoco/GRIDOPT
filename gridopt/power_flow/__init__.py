#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from nr_pf import NRPF
from dc_pf import DCPF
from dc_opf import DCOPF
from augl_pf import AugLPF
from augl_opf import AugLOPF
from method_error import PFmethodError

def new_method(name):

    # Power flow
    if name == 'NRPF':
        return NRPF()
    elif name == 'DCPF':
        return DCPF()
    elif name == 'AugLPF':
        return AugLPF()

    # Optimal power flow
    elif name == 'DCOPF':
        return DCOPF()
    elif name == 'AugLOPF':
        return AugLOPF()

    # invalid
    else:
        raise ValueError('invalid PF method name')
        
    
    
