#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import numpy as np
from types import MethodType
from optalg.opt_solver import OptProblem
from scipy.sparse import eye, bmat, coo_matrix

INF = 1e8

class ProblemTransformer:
    """
    Problem transformer class.
    """

    def __init__(self,p):
        """
        Parameters
        ----------
        p : PFNET Problem
        """
        
        p.analyze()

        self.nx = p.get_num_primal_variables()
        self.m = p.get_num_linear_equality_constraints()
        self.nz = p.G.shape[0]
        self.problem_pfnet = p

    def transform(self):
        """
        Transforms PFNET problem to OPTALG problem.

        Returns
        -------
        new_p : OptProblem
        """

        p = self.problem_pfnet
        nx = self.nx
        nz = self.nz

        new_p = OptProblem()

        new_p.x = np.hstack((p.x.copy(),np.zeros(nz)))
        
        new_p.gphi = np.hstack((p.gphi,np.zeros(nz)))
        new_p.Hphi = coo_matrix((p.Hphi.data,(p.Hphi.row,p.Hphi.col)),shape=(nx+nz,nx+nz))
        
        new_p.A = bmat([[p.A,None],[p.G,-eye(nz)]],format='coo')
        new_p.b = np.hstack((p.b,np.zeros(nz)))
        
        new_p.f = p.f.copy()
        new_p.J = coo_matrix((p.J.data,(p.J.row,p.J.col)),shape=(p.J.shape[0],nx+nz))
        new_p.H_combined = coo_matrix((p.H_combined.data,(p.H_combined.row,p.H_combined.col)),shape=(nx+nz,nx+nz))
        
        new_p.l = np.hstack((-INF*np.ones(nx),p.l))
        new_p.u = np.hstack((INF*np.ones(nx),p.u))
        
        def eval(cls,xz):
            x = xz[:nx]
            z = xz[nx:]
            p.eval(x)
            cls.phi = p.phi
            cls.gphi = np.hstack((p.gphi,np.zeros(nz)))
            cls.Hphi = coo_matrix((p.Hphi.data,(p.Hphi.row,p.Hphi.col)),shape=(nx+nz,nx+nz))
            cls.f = p.f.copy()
            cls.J = coo_matrix((p.J.data,(p.J.row,p.J.col)),shape=(p.J.shape[0],nx+nz))

        def combine_H(cls, coeff, ensure_psd=False):
            p.combine_H(coeff,ensure_psd=ensure_psd)
            cls.H_combined = coo_matrix((p.H_combined.data,(p.H_combined.row,p.H_combined.col)),shape=(nx+nz,nx+nz))
            
        def get_num_primal_variables(cls):
            return nx+nz

        def get_num_linear_equality_constraints(cls):
            return p.get_num_linear_equality_constraints()+nz

        def get_num_nonlinear_equality_constraints(cls):
            return p.get_num_nonlinear_equality_constraints()
    
        new_p.eval = MethodType(eval,new_p)
        new_p.combine_H = MethodType(combine_H,new_p)
        new_p.get_num_primal_variables = MethodType(get_num_primal_variables,new_p)
        new_p.get_num_linear_equality_constraints = MethodType(get_num_linear_equality_constraints,new_p)
        new_p.get_num_nonlinear_equality_constraints = MethodType(get_num_nonlinear_equality_constraints,new_p)
        
        return new_p

    def recover_primal(self,x):

        return x[:self.nx]
    
    def recover_duals(self,lam,nu,mu,pi):
        
        return lam[:self.m],nu,mu[self.nx:],pi[self.nx:]
        

