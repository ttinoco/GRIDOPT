#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from .nr_pf import NRPF
from .dc_pf import DCPF
from .dc_opf import DCOPF
from .augl_pf import AugLPF
from .augl_opf import AugLOPF
from .method import PFmethod
from .method_error import PFmethodError

methods = [NRPF,DCPF,DCOPF,AugLPF,AugLOPF]

def new_method(name):
    """
    Creates a power flow or optimal power flow
    method.
    
    Parameters
    ----------
    name : string
    """
    
    try:
        return methods[list([x.name for x in methods]).index(name)]()
    except ValueError:
        raise ValueError('invalid PF method name')
