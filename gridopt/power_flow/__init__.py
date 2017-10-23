#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from .dc_pf import DCPF
from .dc_opf import DCOPF
from .ac_pf import ACPF
from .ac_opf import ACOPF
from .method import PFmethod
from .method_error import PFmethodError

methods = [DCPF,DCOPF,ACPF,ACOPF]

def new_method(name):
    """
    Creates a power flow or optimal power flow method.
    
    Parameters
    ----------
    name : {``'DCPF'``, ``'DCOPF'``, ``'ACPF'``, ``'ACOPF'``}
    """
    
    try:
        return methods[list([x.name for x in methods]).index(name)]()
    except ValueError:
        raise ValueError('invalid method name')
