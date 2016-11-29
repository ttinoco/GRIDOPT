#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np

def ApplyFunc(args):

    cls = args[0]
    fnc = args[1]
    args = args[2:]
    return getattr(cls,fnc)(*args)

