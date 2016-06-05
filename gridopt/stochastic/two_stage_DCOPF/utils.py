#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np

def ApplyFunc(args):

    cls = args[0]
    fnc = args[1]
    args = args[2:]
    return getattr(cls,fnc)(*args)

class Cache:

    def __init__(self,problem,num_clusters,rate):
        
        self.rate = rate
        self.num_p = problem.num_p
        self.num_w = problem.num_w
        self.num_r = problem.num_r
        self.num_eq = problem.num_eq
        self.problem = problem
        self.num_clusters = num_clusters
        
        self.centroids = np.zeros((0,problem.num_r))
        for i in range(num_clusters):
            p = np.zeros(problem.num_p)
            r = problem.sample_w()
            self.centroids = np.vstack((self.centroids,
                                        r))
        
        self.num_x = self.num_p+self.num_w+self.num_r 
        self.prim_duals = {}
        for i in range(num_clusters):
            x = np.zeros(self.num_x)
            d = np.zeros(self.num_eq+2*self.num_x)
            self.prim_duals[i] = np.hstack((x,d))
        
        assert(0 < rate < 1.)

    def cluster(self,p,r):

        pr = r
        dist = np.square(self.centroids-pr).sum(axis=1)
        i = np.argmin(dist)
        self.centroids[i,:] += self.rate*(pr-self.centroids[i,:])
        xd = self.prim_duals[i]
        return i,xd[:self.num_x],xd[self.num_x:]
        
    def update(self,i,x,d):

        xd = np.hstack((x,d))
        print(np.linalg.norm(xd-self.prim_duals[i]))
        self.prim_duals[i] += self.rate*(xd-self.prim_duals[i])
