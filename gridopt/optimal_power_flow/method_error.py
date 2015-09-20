#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

class OPFmethodError(Exception):
    
    def __init__(self,method,msg):
        method.set_status('error')
        method.set_error_msg(msg)
        self.value = msg
            
    def __str__(self):
        return str(self.value)

class OPFmethodError_BadProblem(OPFmethodError):    
    def __init__(self,method):
        OPFmethodError.__init__(self,method,'error while creating OPF problem')

class OPFmethodError_BadParam(OPFmethodError):
    def __init__(self,method,param=''):
        OPFmethodError.__init__(self,method,'invalid method parameter %s' %param)

class OPFmethodError_ParamNotBool(OPFmethodError):
    def __init__(self,method):
        OPFmethodError.__init__(self,method,'parameter value must be True or False')

class OPFmethodError_SolverError(OPFmethodError):
    def __init__(self,method,msg):
        OPFmethodError.__init__(self,method,msg)
