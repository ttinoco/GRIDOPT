#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

class TS_DCOPF_Method:

    def __init__(self):
        """
        Two-stage DC OPF method class.
        """
        
        #: Results dictionary
        self.results = {}

        #: Parameters dictionary
        self.parameters = {}

        #: Problem
        self.problem = None

    def create_problem(self,net):
        """
        Creates optimization problem.

        Parameters
        ----------
        net : PFNET Network

        Returns
        -------
        prob : OPTALG StochObj_Problem
        """
        
        return None

    def set_parameters(self,parameters):
        """
        Sets solver parameters.

        Parameters
        ----------
        parameters : dict
        """
        
        for key,value in list(parameters.items()):
            if key in self.parameters:
                self.parameters[key] = value

    def get_problem(self):
        """
        Gets problem.
        
        Returns
        -------
        problem
        """
        
        return self.problem

    def get_results(self):
        """
        Gets dictionary with results.

        Returns
        -------
        results : dict
        """

        return self.results

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
