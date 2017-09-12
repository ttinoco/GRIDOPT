#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from .method_error import *

class PFmethod:

    def __init__(self):
        """
        Power flow method class.
        """
        
        self._parameters = {}
        
        self.results = {'status': 'unknown',              # solver status
                        'error msg': '',                  # solver error message
                        'iterations': 0,                  # solver iterations
                        'primal variables': np.zeros(0),  # primal variables
                        'dual variables': 4*[None],       # dual variables
                        'net properties': {},             # network properties
                        'problem': None}                  # PFNET problem

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
        
        self.results['error msg'] = msg

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
        
        self.results['primal variables'] = x

    def set_dual_variables(self,d):
        """
        Sets dual variables.

        Parameters
        ----------
        d : list
        """
        
        self.results['dual variables'] = d

    def set_net_properties(self,np):
        """
        Sets network properties.

        Parameters
        ----------
        np : dictionary
        """
        
        self.results['net properties'] = np

    def set_problem(self,p):
        """
        Sets problem.

        Parameters
        ----------
        p : PFNET problem
        """
        
        self.results['problem'] = p

    def get_parameters(self):
        """
        Gets method parameters.

        Returns
        -------
        params : dict
        """

        return self._parameters

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

        OPTSOLVER_PARAMS = 'optsolver_parameters'

        # Method and solver parameters
        dict_list = [self._parameters]
        if OPTSOLVER_PARAMS in self._parameters:
            dict_list += list(self._parameters[OPTSOLVER_PARAMS].values())
            
        # Parameters
        if params:

            for key,value in list(params.items()):
                if key == OPTSOLVER_PARAMS:
                    continue
                valid_key = False
                for parameter_dict in dict_list:
                    if key in parameter_dict:
                        valid_key = True
                        parameter_dict[key] = value
                if not valid_key:
                    raise PFmethodError_BadParam(self,param=key)

            if OPTSOLVER_PARAMS in params and OPTSOLVER_PARAMS in self._parameters:
                optsolver_params = params[OPTSOLVER_PARAMS]
                for solver_name in self._parameters[OPTSOLVER_PARAMS].keys():
                    if solver_name in optsolver_params:
                        self._parameters[OPTSOLVER_PARAMS][solver_name].update(optsolver_params[solver_name])
                
        # String-based parameters (from command-line utility)
        if strparams:
            for key,valuestr in list(strparams.items()):
                valid_key = False
                for parameter_dict in dict_list:
                    if key in parameter_dict:
                        valid_key = True
                        value = parameter_dict[key]
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
                        else:
                            new_value = valuestr
                        parameter_dict[key] = new_value
                if not valid_key:
                    raise PFmethodError_BadParam(self,param=key)

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
