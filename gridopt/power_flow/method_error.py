#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

class PFmethodError(Exception):
    
    def __init__(self,method,msg):
        method.set_status('error')
        method.set_error_msg(msg)
        self.value = msg
            
    def __str__(self):
        return str(self.value)

class PFmethodError_BadProblem(PFmethodError):    
    def __init__(self,method):
        PFmethodError.__init__(self,method,'error while creating PF problem')

class PFmethodError_BadParam(PFmethodError):
    def __init__(self,method,param=''):
        PFmethodError.__init__(self,method,'invalid method parameter %s' %param)

class PFmethodError_ParamNotBool(PFmethodError):
    def __init__(self,method):
        PFmethodError.__init__(self,method,'parameter value must be True or False')

class PFmethodError_SolverError(PFmethodError):
    def __init__(self,method,msg):
        PFmethodError.__init__(self,method,msg)
