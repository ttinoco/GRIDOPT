#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

class PFmethodError(Exception):
    pass

class PFmethodError_NoProblem(PFmethodError):    
    def __init__(self):
        PFmethodError.__init__(self, 'no problem solved')

class PFmethodError_BadProblem(PFmethodError):    
    def __init__(self):
        PFmethodError.__init__(self, 'error while creating PF problem')

class PFmethodError_BadFlowLimits(PFmethodError):
    def __init__(self):
        PFmethodError.__init__(self, 'invalid flow limits')

class PFmethodError_BadVarLimits(PFmethodError):
    def __init__(self):
        PFmethodError.__init__(self, 'invalid variable limits')

class PFmethodError_BadParams(PFmethodError):
    def __init__(self, keys):
        PFmethodError.__init__(self, 'invalid method parameters: {}'.format(keys))

class PFmethodError_BadParamValue(PFmethodError):
    def __init__(self, key):
        PFmethodError.__init__(self, 'invalid value for parameter {}'.format(key))

class PFmethodError_BadOptSolver(PFmethodError):
    def __init__(self, param=''):
        PFmethodError.__init__(self, 'invalid optimization solver %s' %param)

class PFmethodError_ParamNotBool(PFmethodError):
    def __init__(self):
        PFmethodError.__init__(self, 'parameter value must be True or False')

class PFmethodError_TranVReg(PFmethodError):
        def __init__(self, msg):
            PFmethodError.__init__(self, 'error in transformer voltage regulation :%s' %msg)

class PFmethodError_ShuntVReg(PFmethodError):
        def __init__(self, msg):
            PFmethodError.__init__(self, 'error in shunt voltage regulation: %s' %msg)

class PFmethodError_SolverError(PFmethodError):
    def __init__(self, msg):
        PFmethodError.__init__(self, msg)


