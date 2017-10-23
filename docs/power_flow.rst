.. include:: defs.hrst

.. _power_flow:

**********
Power Flow
**********

The Power Flow (PF) problem consists of determining steady-state voltage magnitudes and angles at every bus of the network as well as any unknown generator powers. On the other hand, the Optimal Power Flow (OPF) problem consists of determining generator powers and network control settings that result in the optimal operation of the network according to some measure, *e.g.*, generation cost. In GRIDOPT, methods for solving PF and OPF problems are represented by objects derived from a :class:`method <gridopt.power_flow.method.PFmethod>` base class.

To solve a PF or OPF problem, one first needs to create an instance of a specific method subclass. This is done using the function :func:`new_method() <gridopt.power_flow.new_method>`, which takes as argument the name of an available PF or OPF method. Available methods are the following: 

* :ref:`dc_pf`
* :ref:`dc_opf`
* :ref:`ac_pf`
* :ref:`ac_opf`

The following code sample creates an instance of the :ref:`ac_pf` method::

  >>> import gridopt

  >>> method = gridopt.power_flow.new_method('ACPF')

Once a method has been instantiated, its parameters can be set using the function :func:`set_parameters() <gridopt.power_flow.method.PFmethod.set_parameters>`. This function takes as argument a dictionary with parameter name-value pairs. Valid parameters include parameters of the method, which are described in the sections below, and parameters of the numerical :class:`solver <optalg.opt_solver.opt_solver.OptSolver>` used by the method. The numerical solvers used by the methods of GRIDOPT belong to the Python package |OPTALG|. The following code sample sets a few parameters of the method created above::

  >>> method.set_parameters({'solver': 'nr', 'quiet': True, 'feastol': 1e-4})

After configuring parameters, a method can be used to solve a problem using the function :func:`solve() <gridopt.power_flow.method.PFmethod.solve>`. This function takes as argument a |network| object, as follows::

  >>> import pfnet

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> method.solve(net)

Information about the execution of the method can be obtained from the :data:`results <gridopt.power_flow.method.PFmethod.results>` attribute of the :class:`method <gridopt.power_flow.method.PFmethod>` object. This dictionary of results includes information such as ``'solver status'``, *e.g.*, ``'solved'`` or ``'error'``, any ``'solver message'``, ``'solver iterations'``, a ``'network snapshot'`` reflecting the solution, and others. The following code sample shows how to extract some results::

  >>> results = method.get_results()

  >>> print results['solver status']
  solved

  >>> print results['solver iterations']
  1

  >>> print results['network snapshot'].bus_v_max
  1.09

If desired, one can update the input |Network| object with the solution found by the method. This can be done with the function :func:`update_network() <gridopt.power_flow.method.PFmethod.update_network>`. This routine not only updates the network quantities treated as variables by the method, but also information about the sensitivity of the optimal objective value with respect to perturbations of the constraints. The following code sample updates the power network with the results obtained by the method and shows the resulting maximum active and reactive bus power mismatches in units of MW and MVAr::

  >>> method.update_network(net)

  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.16e-04 5.67e-03
    
.. _dc_pf: 

DCPF
====

This method is represented by an object of type :class:`DCPF <gridopt.power_flow.dc_pf.DCPF>` and solves a DC power flow problem, which is just a linear system of equations representing |ConstraintDCPF| constraints. The system is solved using one of the |LinSolvers| available in |OPTALG|.

.. _dc_opf: 

DCOPF
=====

This method is represented by an object of type :class:`DCOPF <gridopt.power_flow.dc_opf.DCOPF>` and solves a DC optimal power flow problem, which is just a quadratic program that considers |FunctionGEN_COST|, |FunctionLOAD_UTIL|, |ConstraintDCPF|, |ConstraintBOUND|, *e.g.*, generator and load limits, and |ConstraintDC_FLOW_LIM|. For solving the problem, this method uses by default the |IQP| solver from |OPTALG|.

The parameters of this method are the following:

=========================== ===================================================== =========
Name                        Description                                           Default  
=========================== ===================================================== =========
``'thermal_limits'``        Flag for considering branch flow limits               ``False``
``'renewable_curtailment'`` Flag for allowing curtailment of renewable generators ``False``
``'solver'``                OPTALG optimization solver ``{'iqp','augl','ipopt'}`` ``'iqp'``
=========================== ===================================================== =========

The following example illustrates how to solve a DC OPF problem and extract the optimal generation cost::

  >>> method = gridopt.power_flow.new_method('DCOPF')

  >>> method.solve(net)

  >>> print method.results['solver status']
  solved

  >>> method.update_network(net)

  >>> # generation cost ($/hour)
  >>> print net.gen_P_cost
  7642.63

The sensitivity of the optimal objective value with respect to the power balance constraints can be easily extracted from the network buses::

  >>> bus = net.get_bus(4)
  >>> print "bus %2d %.2e" %(bus.index,bus.sens_P_balance)
  bus 4 3.90e+03
  
Similarly, the sensitivity with respect to branch flow limits can be easily extracted from the network branches::

  >>> branch = net.get_branch(6)
  >>> print "branch %2d %.2e %.2e" %(branch.index,
  ...                                branch.sens_P_u_bound,
  ...                                branch.sens_P_l_bound)
  branch 6 0.00e-00 0.00e-00

Lastly, the sensitivity with respect to generator active power limits can be easily extracted from the network generators::

  >>> gen = net.get_generator(2)
  >>> print "gen %2d %.2e %.2e" %(gen.index,
  ...                             gen.sens_P_u_bound,
  ...                             gen.sens_P_l_bound)
  gen  2 6.36e-04 9.87e+01

As the examples show, GRIDOPT and |PFNET| take care of all the details and allow one to extract solution information easily and intuitively from the |Network| components.

.. _ac_pf:

ACPF
====

This method is represented by an object of type :class:`ACPF <gridopt.power_flow.ac_pf.ACPF>` and solves an AC power flow problem. For doing this, it can use the |NR| solver from |OPTALG| together with "switching" heuristics for modeling local controls. Alternatively, it can formulate the problem as an optimization problem with a convex objective function and *complementarity constraints*, *e.g.*, |ConstraintREG_GEN|, |ConstraintREG_TRAN|, and |ConstraintREG_SHUNT|, for modeling local controls, and solve it using the |AUGL|, |INLP|, or |IPOPT| solver available through |OPTALG|. For now, the parameters of this power flow method are the following:

================= ============================================================ ===========
Name              Description                                                  Default  
================= ============================================================ ===========
``'weight_vmag'`` Weight for bus voltage magnitude regularization              ``1e0``
``'weight_vang'`` Weight for bus voltage angle regularization                  ``1e0``
``'weight_pq'``   Weight for generator power regularization                    ``1e-3``
``'weight_t'``    Weight for transformer tap ratio regularization              ``1e-3``
``'weight_b'``    Weight for shunt susceptance regularization                  ``1e-3``
``'limit_gens'``  Flag for enforcing generator reactive power limits           ``True``
``'lock_taps'``   Flag for locking transformer tap ratios                      ``True``
``'lock_shunts'`` Flag for locking swtiched shunts                             ``True``
``'tap_step'``    Tap ratio acceleration factor (NR heuristics)                ``0.5``
``'shunt_step'``  Susceptance acceleration factor (NR heuristics)              ``0.5``
``'dtap'``        Tap ratio perturbation (NR heuristics)                       ``1e-5``
``'dsus'``        Susceptance perturbation (NR heuristics)                     ``1e-5``
``'vmin_thresh'`` Low-voltage threshold                                        ``1e-1``
``'solver'``      OPTALG optimization solver ``{'nr','inlp','augl','ipopt'}``  ``'augl'``
================= ============================================================ ===========

.. _ac_opf: 

ACOPF
=====

This method is represented by an object of type :class:`ACOPF <gridopt.power_flow.ac_opf.ACOPF>` and solves an AC optimal power flow problem. For doing this, it uses the |AUGL|, |INLP|, or |IPOPT| solver from |OPTALG|. By default, it minimizes |FunctionGEN_COST| subject to voltage magnitude limits, generator power limits, *e.g.*, |ConstraintBOUND|, and |ConstraintACPF|. For now, the parameters of this optimal power flow method are the following:

==================== ============================================================ ===========
Name                 Description                                                  Default  
==================== ============================================================ ===========
``'weight_cost'``    Weight for active power generation cost                      ``1e0`` 
``'weight_vmag'``    Weight for bus voltage magnitude regularization              ``0.``
``'weight_vang'``    Weight for bus voltage angle regularization                  ``0.``
``'weight_pq'``      Weight for generator power regularization                    ``0``
``'weight_t'``       Weight for transformer tap ratio regularization              ``0``
``'weight_b'``       Weight for shunt susceptance regularization                  ``0``
``'thermal_limits'`` Flag for considering |ConstraintAC_FLOW_LIM|                 ``False``
``'vmin_thresh'``    Low-voltage threshold                                        ``1e-1``
``'solver'``         OPTALG optimization solver ``{'augl','inlp','ipopt'}``       ``'augl'``
==================== ============================================================ ===========
