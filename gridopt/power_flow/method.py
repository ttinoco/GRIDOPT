#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from method_error import *

class PFmethod:

    def __init__(self):
        """
        Power flow method class.
        """
        
        #: Results dictionary
        self.results = {'status': 'unknown',
                        'error_msg': '',
                        'variables': np.zeros(0),
                        'iterations': 0}

        #: Parameters dictionary
        self.parameters = {}

    def create_problem(self,net):
        """
        Creates optimization problem.

        Parameters
        ----------
        net : PFNET Network

        Returns
        -------
        prob : PFNET Problem
        """
        
        return None

    def get_info_printer(self):
        """
        Gets function for printing information
        about method progress.

        Returns
        -------
        printer : Function
        """

        return lambda solver,header: None

    def get_results(self):
        """
        Gets dictionary with results.

        Returns
        -------
        results : dict
        """

        return self.results

    def set_error_msg(self,msg):
        """
        Sets method error message.

        Parameters
        ----------
        msg : string
        """
        
        self.results['error_msg'] = msg

    def set_status(self,status):
        """
        Sets method status.

        Parameters
        ----------
        status : string
        """
        
        self.results['status'] = status

    def set_parameters(self,params=None,strparams=None):
        """
        Sets method parameters.

        Parameters
        ----------
        params : dict
                 Name-value pairs
        strparams: dict
                   Name-value pairs where value is a string
        """

        if params:
            for key,value in params.items():
                if self.parameters.has_key(key):
                    self.parameters[key] = value
                else:
                    raise PFmethodError_BadParam(self,param=key)
        if strparams:
            for key,valuestr in strparams.items():
                if self.parameters.has_key(key):
                    value = self.parameters[key]
                    if type(value) is float:
                        new_value = float(valuestr)
                    elif type(value) is int:
                        new_value = int(valuestr)
                    elif type(value) is bool:
                        if valuestr == 'True':
                            new_value = True
                        elif valuestr == 'False':
                            new_value = False
                        else:
                            raise PFmethodError_ParamNotBool(self)
                    elif type(value) is str:
                        new_value = str(valuestr)
                    else:
                        raise PFmethodError_BadParam(self)
                    self.parameters[key] = new_value
                else:
                    raise PFmethodError_BadParam(self)

    def set_results(self,results):
        """
        Sets method results.

        Parameters
        ----------
        results : dict
        """

        self.results = results
                
    def solve(self,net):
        """
        Solves power flow problem.
        
        Parameters
        ----------
        net : PFNET Network
        """        
        
        pass
