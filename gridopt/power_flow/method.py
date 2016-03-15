#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
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
        
        #: Results (dictionary)
        self.results = {'status': 'unknown',              # solver status
                        'error_msg': '',                  # solver error message
                        'iterations': 0,                  # solver iterations
                        'primal_variables': np.zeros(0),  # primal variables
                        'dual_variables': 4*[None],       # dual variables
                        'net_properties': {},             # network properties
                        'problem': None}                  # PFNET problem

        #: Parameters (dictionary)
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

    def set_status(self,status):
        """
        Sets method status.

        Parameters
        ----------
        status : string
        """
        
        self.results['status'] = status

    def set_error_msg(self,msg):
        """
        Sets method error message.

        Parameters
        ----------
        msg : string
        """
        
        self.results['error_msg'] = msg

    def set_iterations(self,k):
        """
        Sets method iterations.

        Parameters
        ----------
        k : int
        """
        
        self.results['iterations'] = k

    def set_primal_variables(self,x):
        """
        Sets primal variables.

        Parameters
        ----------
        x : vector
        """
        
        self.results['primal_variables'] = x

    def set_dual_variables(self,d):
        """
        Sets dual variables.

        Parameters
        ----------
        d : list
        """
        
        self.results['dual_variables'] = d

    def set_net_properties(self,np):
        """
        Sets network properties.

        Parameters
        ----------
        np : dictionary
        """
        
        self.results['net_properties'] = np

    def set_problem(self,p):
        """
        Sets problem.

        Parameters
        ----------
        p : PFNET problem
        """
        
        self.results['problem'] = p

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

    def update_network(self,net):
        """
        Updates network with results.

        Parameters
        ----------
        net : PFNET Network
        """

        pass
