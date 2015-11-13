#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import pfnet
import numpy as np
from method_error import *
from method import OPFmethod
from optalg.opt_solver import OptSolverError,OptCallback,OptTermination,OptSolverAugL

class AugLOPF(OPFmethod):
    """
    Augmented Lagrangian-based optimal power flow.
    """
 
    parameters = {'alpha_Pc':1e-2,     # for generation cost
                  'alpha_Pl':1e-4,     # for soft limits
                  'alpha_Pr':1e-5,     # for regularization
                  'vmin_thresh':0.1}   # threshold for vmin
                   
    def __init__(self):

        OPFmethod.__init__(self)
        self.parameters = AugLOPF.parameters.copy()
        self.parameters.update(OptSolverAugL.parameters)

    def create_problem(self,net):

        # Parameters
        params = self.parameters
        alpha_Pc = params['alpha_Pc']
        alpha_Pl = params['alpha_Pl']
        alpha_Pr = params['alpha_Pr']
        
        # Clear flags
        net.clear_flags()
        
        # Voltage magnitudes
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_ANY,
                      pfnet.BUS_VAR_VMAG)
        
        # Voltage angles
        net.set_flags(pfnet.OBJ_BUS,
                      pfnet.FLAG_VARS,
                      pfnet.BUS_PROP_NOT_SLACK,
                      pfnet.BUS_VAR_VANG)

        # Generator active power
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS|pfnet.FLAG_BOUNDED,
                      pfnet.GEN_PROP_ANY,
                      pfnet.GEN_VAR_P)

        # Generator reactive power
        net.set_flags(pfnet.OBJ_GEN,
                      pfnet.FLAG_VARS|pfnet.FLAG_BOUNDED,
                      pfnet.GEN_PROP_REG,
                      pfnet.GEN_VAR_Q)

        try:
            assert(net.num_vars == (2*net.num_buses-net.get_num_slack_buses() +
                                    net.num_gens + 
                                    net.get_num_reg_gens()))
            assert(net.num_bounded == net.num_gens + net.get_num_reg_gens())
        except AssertionError:
            raise OPFmethodError_BadProblem(self)
                                    
        # Set up problem
        problem = pfnet.Problem()
        problem.set_network(net)
        problem.add_constraint(pfnet.CONSTR_TYPE_PF)
        problem.add_constraint(pfnet.CONSTR_TYPE_BOUND) 
        problem.add_function(pfnet.FUNC_TYPE_GEN_COST,alpha_Pc)
        problem.add_function(pfnet.FUNC_TYPE_SLIM_VMAG,alpha_Pl)
        problem.add_function(pfnet.FUNC_TYPE_REG_VANG,alpha_Pr)
        problem.add_function(pfnet.FUNC_TYPE_REG_PQ,alpha_Pr)
        problem.analyze()
        
        # Return
        return problem

    def get_info_printer(self):

        def info_printer(solver,header):
            net = solver.problem.network
            if header:
                print '{0:^5}'.format('vmax'),
                print '{0:^5}'.format('vmin'),
                print '{0:^6}'.format('bvvio'),
                print '{0:^6}'.format('gQvio'),
                print '{0:^6}'.format('gPvio')
            else:
                print '{0:^5.2f}'.format(net.bus_v_max),
                print '{0:^5.2f}'.format(net.bus_v_min),
                print '{0:^6.0e}'.format(net.bus_v_vio),
                print '{0:^6.0e}'.format(net.gen_Q_vio),
                print '{0:^6.0e}'.format(net.gen_P_vio)
        return info_printer
            
    def solve(self,net):
        
        # Parameters
        params = self.parameters
        vmin_thresh = params['vmin_thresh']

        # Problem
        problem = self.create_problem(net)

        # Callbacks

        # Termination
        def t1(s):
            if s.problem.network.bus_v_min < vmin_thresh:
                return True
            else:
                return False
        
        # Info printer
        info_printer = self.get_info_printer()
        
        # Set up solver
        solver = OptSolverAugL()
        solver.set_parameters(params)
        solver.add_termination(OptTermination(t1,'low voltage'))
        solver.set_info_printer(info_printer)
        
        # Solve
        try:
            solver.solve(problem)
        except OptSolverError,e:
            raise OPFmethodError_SolverError(self,e)
        finally:
            
            # Get results
            self.results = {'status': solver.get_status(),
                            'error_msg': solver.get_error_msg(),
                            'variables': solver.get_primal_variables(),
                            'iterations': solver.get_iterations()}
            self.results.update(net.get_properties())
        
            # Show
            problem.show()
