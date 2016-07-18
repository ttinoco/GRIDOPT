#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

class MS_DCOPF_Method:

    def __init__(self):
        """
        Multi-stage DC OPF method class.
        """
        
        #: Results dictionary
        self.results = {}

        #: Parameters dictionary
        self.parameters = {}

    def create_problem(self,net,forecast):
        """
        Creates optimization problem.

        Parameters
        ----------
        net : PFNET Network
        forecast : dict

        Returns
        -------
        prob : OPTALG StochObjMS_Problem
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
                
    def solve(self,net,forecast):
        """
        Solves power flow problem.
        
        Parameters
        ----------
        net : PFNET Network
        forecast : dict
        
        Returns
        -------
        policy : 
        """        
        
        pass
