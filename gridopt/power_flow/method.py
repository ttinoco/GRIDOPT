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
        
        self.results = {'solver name': 'unknown',
                        'solver status': 'unknown',
                        'solver message': '',
                        'solver iterations': 0,
                        'solver time': np.nan,
                        'solver primal variables': None,
                        'solver dual variables': None,
                        'problem' : None,
                        'problem time' : np.nan,
                        'network snapshot' : None}

    def create_problem(self,net):
        """
        Creates optimization problem.

        Parameters
        ----------
        net : |Network|

        Returns
        -------
        prob : |Problem|
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

    def set_solver_name(self, name):
        """
        Sets solver name.

        Parameters
        ----------
        name : string
        """
        
        self.results['solver name'] = name

    def set_solver_status(self, status):
        """
        Sets solver status.

        Parameters
        ----------
        status : string
        """
        
        self.results['solver status'] = status

    def set_solver_message(self, msg):
        """
        Sets solver message.

        Parameters
        ----------
        msg : string
        """
        
        self.results['solver message'] = msg

    def set_solver_iterations(self, k):
        """
        Sets solver iterations.

        Parameters
        ----------
        k : int
        """
        
        self.results['solver iterations'] = k

    def set_solver_time(self, t):
        """
        Sets solver time in seconds.

        Parameters
        ----------
        t : float
        """
        
        self.results['solver time'] = t

    def set_solver_primal_variables(self, x):
        """
        Sets solver primal variables.

        Parameters
        ----------
        x : vector
        """
        
        self.results['solver primal variables'] = x

    def set_solver_dual_variables(self, d):
        """
        Sets solver dual variables.

        Parameters
        ----------
        d : list
        """
        
        self.results['solver dual variables'] = d

    def set_problem(self, p):
        """
        Sets problem.
         
        Parameters
        ----------
        p : |Problem|
        """

        self.results['problem'] = p

    def set_problem_time(self, t):
        """
        Sets problem construction time in seconds.

        Parameters
        ----------
        t : float
        """
        
        self.results['problem time'] = t

    def set_network_snapshot(self, net):
        """
        Sets network snapshot.

        Parameters
        ----------
        net : |Network|
        """
        
        self.results['network snapshot'] = net

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

        invalid_params = []

        SOLVER_PARAMS = 'solver_parameters'

        # List of method/solver parameter dictionaries
        dict_list = [self._parameters]
        if SOLVER_PARAMS in self._parameters:
            dict_list += list(self._parameters[SOLVER_PARAMS].values())
            
        # Parameters
        if params:
            for key,value in list(params.items()):
                if key == SOLVER_PARAMS:
                    continue
                valid_key = False
                for parameter_dict in dict_list:
                    if key in parameter_dict:
                        valid_key = True
                        parameter_dict[key] = value
                if not valid_key:
                    invalid_params.append(key)
                    
            if SOLVER_PARAMS in params and SOLVER_PARAMS in self._parameters:
                solver_params = params[SOLVER_PARAMS]
                for solver_name in self._parameters[SOLVER_PARAMS].keys():
                    if solver_name in solver_params:
                        self._parameters[SOLVER_PARAMS][solver_name].update(solver_params[solver_name])
                
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
                                raise PFmethodError_ParamNotBool()
                        else:
                            new_value = valuestr
                        parameter_dict[key] = new_value
                if not valid_key:
                    invalid_params.append(key)

        # Invalid params
        if invalid_params:
            raise PFmethodError_BadParams(invalid_params)
                
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
        net : |Network|
        """        
        
        pass

    def update_network(self,net):
        """
        Updates network with results.

        Parameters
        ----------
        net : |Network|
        """

        if self.results['network snapshot'] is not None:
            net.copy_from_network(self.results['network snapshot'])
            net.update_properties()

