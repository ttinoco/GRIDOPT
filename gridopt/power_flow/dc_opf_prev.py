#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2016, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from method_error import *
from method import PFmethod
from scipy.sparse import triu,coo_matrix,bmat,eye
from optalg.opt_solver import OptSolverError,OptSolverIQP,QuadProblem

class DCOPF_Prev(PFmethod):
    """
    Preventive DC optimal power flow method.
    """

    name = 'DCOPF_Prev'

    parameters = {'quiet' : False,
                  'thermal_limits': True,
                  'thermal_factor': 1.,
                  'inf_flow': 1e4}
    
    def __init__(self):
        
        PFmethod.__init__(self)
        self.parameters = DCOPF_Prev.parameters.copy()
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

        try:
            num_gvar =  len([g for g in net.generators if 
                             (not g.is_on_outage()) and g.is_P_adjustable()])
            assert(net.num_bounded == num_gvar)
            assert(net.num_vars == net.num_buses-net.get_num_slack_buses()+num_gvar)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
            
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_LBOUND)
        problem.add_constraint(pfnet.CONSTR_TYPE_DCPF)
        problem.add_constraint(pfnet.CONSTR_TYPE_DC_FLOW_LIM)
        problem.add_function(pfnet.FUNC_TYPE_GEN_COST,1.)
        problem.analyze()
        
        # Return
        return problem
            
    def solve(self,net,contingencies):
        
        # Parameters
        params = self.parameters
        thermal_limits = params['thermal_limits']
        thermal_factor = params['thermal_factor']
        inf_flow = params['inf_flow']

        # Problem (base)
        problem = self.create_problem(net)
        
        # Projections
        Pw = net.get_var_projection(pfnet.OBJ_BUS,pfnet.BUS_VAR_VANG)
        Pp = net.get_var_projection(pfnet.OBJ_GEN,pfnet.GEN_VAR_P)
      
        # Construct QP
        x = problem.get_init_point()
        p = Pp*x
        w = Pw*x

        c_flows = problem.find_constraint(pfnet.CONSTR_TYPE_DC_FLOW_LIM)
        c_bounds = problem.find_constraint(pfnet.CONSTR_TYPE_LBOUND)
       
        problem.eval(x)

        phi = problem.phi
        
        Hp = Pp*(problem.Hphi + problem.Hphi.T - triu(problem.Hphi))*Pp.T
        gp = Pp*problem.gphi - Hp*p
        
        G = problem.A*Pp.T
        W = -problem.A*Pw.T
        b = problem.b.copy()
        
        lz = c_flows.l.copy()
        uz = c_flows.u.copy()
        J = c_flows.G*Pw.T

        # Flow limit expansion
        dz = (thermal_factor-1.)*(uz-lz)/2.
        if not thermal_limits:
            dz += inf_flow
        lz -= dz
        uz += dz

        lw = Pw*c_bounds.l
        uw = Pw*c_bounds.u
        Iw = Pw*c_bounds.G*Pw.T

        lp = Pp*c_bounds.l
        up = Pp*c_bounds.u
        Ip = Pp*c_bounds.G*Pp.T

        GWJ_list = [(G,W,J)]
        u_list = [up,uw,uz]
        l_list = [lp,lw,lz]
        b_list = [b,np.zeros(J.shape[0])]
        nz_list = [J.shape[0]]
        
        for cont in contingencies:

            # apply contingency
            cont.apply()

            problem.analyze()

            G = problem.A*Pp.T
            W = -problem.A*Pw.T
            b = problem.b.copy()
            
            lz = c_flows.l.copy()
            uz = c_flows.u.copy()
            J = c_flows.G*Pw.T

            # Flow limit expansion
            dz = (thermal_factor-1.)*(uz-lz)/2.
            if not thermal_limits:
                dz += inf_flow
            lz -= dz
            uz += dz
            
            GWJ_list.append((G,W,J))
            u_list += [uw,uz]
            l_list += [lw,lz]
            b_list += [b,np.zeros(J.shape[0])]
            nz_list.append(J.shape[0])

            # clear contingency
            cont.clear()
          
        problem.analyze()
  
        A = []
        num_blocks = len(GWJ_list)
        for i in range(num_blocks):

            G,W,J = GWJ_list[i]

            row1 = (2*num_blocks+1)*[None]
            row1[0] = G
            row1[2*i+1] = -W
            A.append(row1)

            row2 = (2*num_blocks+1)*[None]
            row2[2*i+1] = J
            row2[2*i+2] = -eye(J.shape[0],format='coo')
            A.append(row2)

        A = bmat(A,format='coo')
        b = np.hstack((b_list))
        l = np.hstack((l_list))
        u = np.hstack((u_list))

        n = A.shape[1]
        ng = Pp.shape[0]
        nw = Pw.shape[0]
        nz = nz_list[0]
        nr = n-ng
        m = A.shape[0]

        Zr = coo_matrix((nr,nr))
        zr = np.zeros(nr)
        
        H = bmat([[Hp,None],[None,Zr]],format='coo')/net.base_power # scaled
        g = np.hstack((gp,zr))/net.base_power                       # scaled

        y = np.hstack((p,zr))

        # Check limits
        if not np.all(l < u):
            raise PFmethodError_BadFlowLimits(self)
        
        # Other checks
        try:
            assert(ng+nw == net.num_vars)
            assert(b.shape == (m,))
            assert((Ip-eye(Pp.shape[0])).nnz == 0)
            assert((Iw-eye(Pw.shape[0])).nnz == 0)
            assert(l.shape == (n,))
            assert(u.shape == (n,))
            assert(np.abs(phi-net.base_power*(0.5*np.dot(y,H*y)+np.dot(g,y))) < 1e-8)
            assert(H.shape == (n,n))
            assert(m == num_blocks*net.num_buses+sum(nz_list))
            assert(np.linalg.norm(x-Pp.T*p-Pw.T*w,np.inf) < 1e-8)
        except AssertionError:
            raise PFmethodError_BadProblem(self)
                        
        QPproblem = QuadProblem(H,g,A,b,l,u)
        
        # Set up solver
        solver = OptSolverIQP()
        solver.set_parameters(params)
        
        # Solve
        try:
            solver.solve(QPproblem)
        except OptSolverError,e:
            raise PFmethodError_SolverError(self,e)
        finally:
            
            # Update net properties
            pwz = solver.get_primal_variables()
            x = Pp.T*pwz[:ng]+Pw.T*pwz[ng:ng+nw]
            z = pwz[ng+nw:ng+nw+nz]
            net.update_properties(x)

            # Prepare duals
            lam,nu,mu,pi = solver.get_dual_variables()
            lam = lam[:net.num_buses+nz]
            mu_p = mu[:ng]
            mu_w = mu[ng:ng+nw]
            mu_z = mu[ng+nw:ng+nw+nz]
            mu = np.hstack((Pp.T*mu_p+Pw.T*mu_w,mu_z))            
            pi_p = pi[:ng]
            pi_w = pi[ng:ng+nw]
            pi_z = pi[ng+nw:ng+nw+nz]
            pi = np.hstack((Pp.T*pi_p+Pw.T*pi_w,pi_z))
            
            # Get results
            self.set_status(solver.get_status())
            self.set_error_msg(solver.get_error_msg())
            self.set_iterations(solver.get_iterations())
            self.set_primal_variables(np.hstack((x,z)))
            self.set_dual_variables([lam,nu,mu,pi])
            self.set_net_properties(net.get_properties())
            self.set_problem(problem)
            
    def update_network(self,net):
    
        # Get data
        problem = self.results['problem']
        xz = self.results['primal_variables']
        lam,nu,mu,pi = self.results['dual_variables']
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
