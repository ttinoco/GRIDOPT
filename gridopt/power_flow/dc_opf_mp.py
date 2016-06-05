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
from scipy.sparse import triu,coo_matrix,bmat,eye,block_diag
from optalg.opt_solver import OptSolverError,OptSolverIQP,QuadProblem

class DCOPF_MP(PFmethod):
    """
    Multi-period DC optimal power flow method.
    """

    name = 'DCOPF_MP'
        
    parameters = {'quiet' : False,
                  'thermal_limits': True,
                  'thermal_factor': 1.,
                  'fixed_total_load': False,
                  'inf_flow': 1e4}
                                    
    def __init__(self):

        PFmethod.__init__(self)
        self.parameters = DCOPF_MP.parameters.copy()
        self.parameters.update(OptSolverIQP.parameters)

    def create_problem(self,net):
        
        # Parameters
        params = self.parameters
        
        # Clear flags
        net.clear_flags()
        
        # Set flags
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VANG)
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS|pfnet.FLAG_BOUNDED,
                      pfnet.GEN_PROP_P_ADJUST|pfnet.GEN_PROP_NOT_OUT,
                      pfnet.GEN_VAR_P)
        net.set_flags(pfnet.OBJ_LOAD,
                      pfnet.FLAG_VARS|pfnet.FLAG_BOUNDED,
                      pfnet.LOAD_PROP_P_ADJUST,
                      pfnet.LOAD_VAR_P)

        try:
            num_gvar =  len([g for g in net.generators if 
                             (not g.is_on_outage()) and g.is_P_adjustable()])
            assert(net.num_bounded == num_gvar + net.get_num_P_adjust_loads())
            assert(net.num_vars == (net.num_buses-net.get_num_slack_buses()+
                                    num_gvar+net.get_num_P_adjust_loads()))
        except AssertionError:
            raise PFmethodError_BadProblem(self)
            
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_LBOUND)
        problem.add_constraint(pfnet.CONSTR_TYPE_DCPF)
        problem.add_constraint(pfnet.CONSTR_TYPE_DC_FLOW_LIM)
        problem.add_function(pfnet.FUNC_TYPE_GEN_COST,1.)
        problem.add_function(pfnet.FUNC_TYPE_LOAD_UTIL,-1.)
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net,T,net_modifier):
        """
        Solves multi-period DC OPF problem.

        Parameters
        ----------
        net : Network
        T : int (time horizon)
        net_modifier : function(net,t)
        """
        
        # Parameters
        params = self.parameters
        thermal_limits = params['thermal_limits']
        thermal_factor = params['thermal_factor']
        fixed_total_load = params['fixed_total_load']
        inf_flow = params['inf_flow']

        # Check parameters
        if thermal_factor < 0.:
            raise PFmethodError_BadParam(self,'thermal_factor')
       
        # Construct MP problem
        data = []
        Lproj = []
        Ltot = 0
        for t in range(T):
            
            # Modify network
            net_modifier(net,t)
            
            # Create problem
            problem = self.create_problem(net)

            # Express in standard form
            x = problem.get_init_point()
            problem.eval(x)
            c_flows = problem.find_constraint(pfnet.CONSTR_TYPE_DC_FLOW_LIM)
            c_bounds = problem.find_constraint(pfnet.CONSTR_TYPE_LBOUND)
                       
            Hx = problem.Hphi + problem.Hphi.T - triu(problem.Hphi)
            gx = problem.gphi - Hx*x
        
            Ax = problem.A
            bx = problem.b
            
            lz = c_flows.l
            uz = c_flows.u
            Gz = c_flows.G
            
            dz = (thermal_factor-1.)*(uz-lz)/2.
            if not thermal_limits:
                dz += inf_flow
            lz -= dz
            uz += dz
  
            lx = c_bounds.l
            ux = c_bounds.u
            Gx = c_bounds.G
            
            nx = net.num_vars
            nz = net.get_num_branches_not_on_outage()
            n = nx+nz
      
            Iz = eye(nz)
            Oz = coo_matrix((nz,nz))
            oz = np.zeros(nz)
       
            if fixed_total_load:
                Pl = net.get_var_projection(pfnet.OBJ_LOAD,pfnet.LOAD_VAR_P)
                Lproj.append(bmat([[Pl,coo_matrix((Pl.shape[0],nz))]]))
                Ltot = Ltot + Pl*x 
 
            H = bmat([[Hx,None],[None,Oz]],format='coo')/net.base_power # scaled
            g = np.hstack((gx,oz))/net.base_power                       # scaled
            
            A = bmat([[Ax,None],[Gz,-Iz]],format='coo')
            b = np.hstack((bx,oz))
            
            l = np.hstack((lx,lz))
            u = np.hstack((ux,uz))
            
            y = np.hstack((x,oz))
            
            if not np.all(lz < uz):
                raise PFmethodError_BadFlowLimits(self)
                
            if not np.all(lx < ux):
                raise PFmethodError_BadVarLimits(self)

            try:
                assert(Gx.shape == (nx,nx))
                assert(np.all(Gx.row == Gx.col))
                assert(np.all(Gx.data == np.ones(nx)))
                assert(Gz.shape == (net.get_num_branches_not_on_outage(),nx))
                assert(l.shape == (n,))
                assert(u.shape == (n,))
                assert(np.all(l < u))
                assert(np.abs(problem.phi-net.base_power*(0.5*np.dot(y,H*y)+np.dot(g,y))) < 1e-7)
                assert(H.shape == (n,n))
                assert(A.shape == (net.num_buses+nz,n))
            except AssertionError:
                raise PFmethodError_BadProblem(self)
           
            data.append((H,g,A,b,l,u))
            
        # Construct QP
        H,g,A,b,l,u = list(zip(*data))
        H = block_diag(H,format='coo')
        g = np.hstack(g)
        A = block_diag(A,format='coo')
        b = np.hstack(b)
        l = np.hstack(l)
        u = np.hstack(u)
        if fixed_total_load:
            L = bmat([Lproj])
            A = bmat([[A],[L]])
            b = np.hstack((b,Ltot))
        QPproblem = QuadProblem(H,g,A,b,l,u)
        
        # Set up solver
        solver = OptSolverIQP()
        solver.set_parameters(params)
        
        # Solve QP
        try:
            solver.solve(QPproblem)
        except OptSolverError as e:
            raise PFmethodError_SolverError(self,e)
        finally:
            
            # Separate results
            offset_n = 0
            offset_A = 0
            duals = []
            primals = []
            problems = []
            properties = []
            x = solver.get_primal_variables()
            lam,nu,mu,pi = solver.get_dual_variables()
            for t in range(T):
                
                # Data
                H,g,A,b,l,u = data[t]

                # Modify net
                net_modifier(net,t)
                
                # Create problem
                problems.append(self.create_problem(net))
                
                # Primals
                primals.append(x[offset_n:offset_n+g.size])

                # Duals
                duals.append((lam[offset_A:offset_A+A.shape[0]],
                              None,
                              mu[offset_n:offset_n+g.size],
                              pi[offset_n:offset_n+g.size]))
        

                # Properties
                net.update_properties(x[offset_n:offset_n+net.num_vars])
                properties.append(net.get_properties())
            
                # Offsets
                offset_n += g.size
                offset_A += A.shape[0]
 
            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(primals)
            self.set_dual_variables(duals)
            self.set_net_properties(properties)
            self.set_problem(problems)
            
    def update_network(self,net,t,net_modifier):
        """
        Updates the network with part of the solution that 
        corresponds to the given time. 

        Parameters
        ----------
        net : Network
        T : int (time)
        net_modifier : function(net,t)
        """
   
        # Modify network
        net_modifier(net,t)
 
        # Get data
        problem = self.create_problem(net)
        xz = self.results['primal_variables'][t]
        lam,nu,mu,pi = self.results['dual_variables'][t]
        nx = net.num_vars
        nz = net.get_num_branches_not_on_outage()
        n = nx+nz
        
        # No problem
        if problem is None:
            raise PFmethodError_NoProblem(self)
 
        # Checks
        assert(problem.x.shape == (nx,))
        assert(problem.A.shape == (net.num_buses,nx))
        assert(problem.G.shape == (nx+nz,nx))
        assert(xz.shape == (nx+nz,))
        assert(lam.shape == (net.num_buses+nz,))
        assert(mu.shape == (n,))
        assert(pi.shape == (n,))
        assert(nu is None)

        # Network quantities
        net.set_var_values(xz[:nx])

        # Network properties
        net.update_properties()
        
        # Network sensitivities
        net.clear_sensitivities()
        problem.store_sensitivities(lam[:net.num_buses]*net.base_power,
                                    nu,
                                    mu*net.base_power,
                                    pi*net.base_power)
