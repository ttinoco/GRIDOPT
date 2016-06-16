#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np

class Node:

    def __init__(self,w,parent=):
        self.w = w
        self.children = []
        self.parent = 

    def get_w(self):
        return self.w

    def get_children(self):
        return self.children

    def add_child(self,child):
        self.children.append(child)

class ScenarioTree:

    def __init__(self,problem,branching_factor,seed=None):
        """
        Creates scenario tree for multistage
        stochastic optimization problem.
        
        Parameters
        ----------
        problem : 
        branching_factor :
        seed :
        """

        if seed is not None:
            np.random.seed(seed)
      
        T = problem.get_num_stages()
  
        root = Node(problem.sample_w(0,[]))
        nodes = [root]
        branch = []
        
        for t in range(1,T):

            new_nodes = []
            for node in nodes:

                branch.append(node)

                for i in range(branching_factor):
                    child = Node(problem.sample_w(t,map(lamnbda n: n.get_w(),branch)))
                    node.add_child(child)
                    
            
            
