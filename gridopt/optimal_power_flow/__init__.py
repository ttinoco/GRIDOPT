#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from dc import DCOPF

def new_method(name):
    
    if name == 'DCOPF':
        return DCOPF()
    else:
        raise ValueError('invalid OPF method name')
        
    
    
