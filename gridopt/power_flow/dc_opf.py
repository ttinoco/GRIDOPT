#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from .method_error import *
from .method import PFmethod
from numpy.linalg import norm
from scipy.sparse import triu,coo_matrix,bmat,eye
from optalg.opt_solver import OptSolverError,OptSolverIQP,QuadProblem

class DCOPF(PFmethod):
    """
    DC optimal power flow method.
    """

    name = 'DCOPF'
        
    parameters = {'quiet' : False,
                  'thermal_limits': True,
                  'thermal_factor': 1.,
                  'inf_flow': 1e4,
                  'vargen_curtailment': False}
                                    
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = DCOPF.parameters.copy()
        self.parameters.update(OptSolverIQP.parameters)

    def create_problem(self,net):
        
        # Parameters
        params = self.parameters
        
        # Clear flags
        net.clear_flags()
        
        # Set flags
        net.set_flags('bus',
                      'variable',
                      'not slack',
                      'voltage angle')
        net.set_flags('generator',
                      ['variable','bounded'],
                      ['adjustable active power','not on outage'],
                      'active power')
        net.set_flags('load',
                      ['variable','bounded'],
                      'adjustable active power',
                      'active power')
        if params['vargen_curtailment']:
            net.set_flags('variable generator',
                          ['variable','bounded'],
                          'any',
                          'active power')

        try:
            num_gvar =  len([g for g in net.generators if 
                             (not g.is_on_outage()) and g.is_P_adjustable()])
            num_cur = net.num_var_generators if params['vargen_curtailment'] else 0
            assert(net.num_bounded == (num_gvar+net.get_num_P_adjust_loads()+num_cur)*net.num_periods)
            assert(net.num_vars == (net.num_buses-net.get_num_slack_buses()+
                                    num_gvar+net.get_num_P_adjust_loads()+
                                    num_cur)*net.num_periods)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
            
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint('variable bounds')
        problem.add_constraint('DC power balance')
        problem.add_constraint('DC branch flow limits')
        problem.add_function('generation cost',1.)
        problem.add_function('consumption utility',-1.)
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        thermal_limits = params['thermal_limits']
        thermal_factor = params['thermal_factor']
        inf_flow = params['inf_flow']

        # Check parameters
        if thermal_factor < 0.:
            raise PFmethodError_BadParam(self,'thermal_factor')

        # Problem
        problem = self.create_problem(net)

        # Renewables
        Pr = net.get_var_projection('variable generator','active power')
       
        # Construct QP
        x = problem.get_init_point()
        problem.eval(x)
        c_flows = problem.find_constraint('DC branch flow limits')
        c_bounds = problem.find_constraint('variable bounds')
        
        Hx = problem.Hphi + problem.Hphi.T - triu(problem.Hphi)
        gx = problem.gphi - Hx*x
        
        Ax = problem.A
        bx = problem.b
        
        lz = c_flows.l
        uz = c_flows.u
        Gz = c_flows.G

        # Flow limit expansion
        dz = (thermal_factor-1.)*(uz-lz)/2.
        if not thermal_limits:
            dz += inf_flow
        lz -= dz
        uz += dz
  
        lx = c_bounds.l
        ux = c_bounds.u
        Gx = c_bounds.G
        
        ux += Pr.T*Pr*(x-ux) # correct limit for curtailment

        nx = net.num_vars
        nz = net.get_num_branches_not_on_outage()*net.num_periods
        n = nx+nz

        Iz = eye(nz)
        Oz = coo_matrix((nz,nz))
        oz = np.zeros(nz)
        
        H = bmat([[Hx,None],[None,Oz]],format='coo')
        g = np.hstack((gx,oz))

        A = bmat([[Ax,None],[Gz,-Iz]],format='coo')
        b = np.hstack((bx,oz))

        l = np.hstack((lx,lz))
        u = np.hstack((ux,uz))

        y = np.hstack((x,oz))

        # Check flow limits
        if not np.all(lz < uz):
            raise PFmethodError_BadFlowLimits(self)
        
        # Check variable limits
        if not np.all(lx < ux):
            raise PFmethodError_BadVarLimits(self)

        # Other checks
        try:
            assert(Gx.shape == (nx,nx))
            assert(np.all(Gx.row == Gx.col))
            assert(np.all(Gx.data == np.ones(nx)))
            assert(Gz.shape == (net.get_num_branches_not_on_outage()*net.num_periods,nx))
            assert(l.shape == (n,))
            assert(u.shape == (n,))
            assert(np.all(l < u))
            assert(norm(np.hstack((problem.gphi,oz))-(H*y+g)) < 1e-10*(1.+norm(problem.gphi)))
            assert(H.shape == (n,n))
            assert(A.shape == (net.num_buses*net.num_periods+nz,n))
        except AssertionError:
            raise PFmethodError_BadProblem(self)
            
        QPproblem = QuadProblem(H,g,A,b,l,u)
        
        # Set up solver
        solver = OptSolverIQP()
        solver.set_parameters(params)
        
        # Solve
        try:
            solver.solve(QPproblem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:

            # Update net properties
            net.update_properties(solver.get_primal_variables()[:net.num_vars])
            
            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(solver.get_primal_variables())
            self.set_dual_variables(solver.get_dual_variables())
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)
            
    def update_network(self,net):
    
        # Get data
        problem = self.results['problem']
        xz = self.results['primal_variables']
        lam,nu,mu,pi = self.results['dual_variables']
        nx = net.num_vars
        nz = net.get_num_branches_not_on_outage()*net.num_periods
        n = nx+nz
        
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == (nx,))
        assert(problem.A.shape == (net.num_buses*net.num_periods,nx))
        assert(problem.G.shape == (nx+nz,nx))
        assert(xz.shape == (nx+nz,))
        assert(lam.shape == (net.num_buses*net.num_periods+nz,))
        assert(mu.shape == (n,))
        assert(pi.shape == (n,))
        assert(nu is None or not nu.size)

        # Network quantities
        net.set_var_values(xz[:nx])

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam[:net.num_buses*net.num_periods],
                                    nu,
                                    mu,
                                    pi)
