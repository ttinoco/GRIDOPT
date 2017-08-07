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

Once a method has been instantiated, its parameters can be set using the function :func:`set_parameters() <gridopt.power_flow.method.PFmethod.set_parameters>`. This function takes as argument a dictionary with parameter name and value pairs. Valid parameters include parameters of the method, which are described in the sections below, and parameters of the numerical solver used by the method. The numerical solvers used by the methods of GRIDOPT belong to the Python package `OPTALG`_. The following code sample sets a few parameters of the method created above::

  >>> method.set_parameters({'optsolver': 'nr', 'quiet': True, 'feastol': 1e-4})

After configuring parameters, a method can be used to solve a problem using the function :func:`solve() <gridopt.power_flow.method.PFmethod.solve>`. This function takes as argument a `network`_ object, as follows::

  >>> import pfnet

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> method.solve(net)

Information about the execution of the method can be obtained from the :data:`results <gridopt.power_flow.method.PFmethod.results>` attribute of the :class:`method <gridopt.power_flow.method.PFmethod>` object. This dictionary of results includes information such as ``'status'``, *e.g.*, ``'solved'`` or ``'error'``, any error message (``'error_msg'``), solver ``'iterations'``, the `optimization problem`_ (``'problem'``) constructed, and `network properties`_ at the point found by the method (``'net properties'``). The following code sample shows how to extract some results::

  >>> results = method.get_results()

  >>> print results['status']
  solved

  >>> print results['iterations']
  1

  >>> problem = results['problem']
  >>> problem.show()
  
  Problem
  functions  : 0
  constraints: 4
      AC power balance
      generator active power participation
      generator reactive power participation
      variable fixing

  >>> print results['net properties']['bus_v_max']
  1.09

If desired, one can update the `network`_ object with the solution found by the method. This can be done with the function :func:`update_network() <gridopt.power_flow.method.PFmethod.update_network>`. This routine not only updates the network quantities treated as variables by the method, but also information about the sensitivity of the optimal objective value with respect to perturbations of the constraints. The following code sample updates the power network with the results obtained by the method and shows the resulting maximum active and reactive bus power mismatches in units of MW and MVAr::

  >>> method.update_network(net)

  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.16e-04 5.67e-03
    
.. _dc_pf: 

DCPF
====

This method is represented by an object of type :class:`DCPF <gridopt.power_flow.dc_pf.DCPF>` and solves a DC power flow problem, which is just a linear system of equations representing `DC power balance constraints`_.  The system is solved using one of the ``linear solvers`` available in `OPTALG`_.

.. _dc_opf: 

DCOPF
=====

This method is represented by an object of type :class:`DCOPF <gridopt.power_flow.dc_opf.DCOPF>` and solves a DC optimal power flow problem, which is just a quadratic program that considers `active power generation cost`_, `active power consumption utility`_, `DC power balance constraints`_, `variable limits`_, *e.g.*, generator and load limits, and `DC branch flow limits`_. For solving the problem, this method uses by default the `IQP`_ solver from `OPTALG`_.

The parameters of this method are the following:

=========================== ===================================================== =========
Name                        Description                                           Default  
=========================== ===================================================== =========
``'thermal_limits'``        Flag for considering branch flow limits               ``False``
``'renewable_curtailment'`` Flag for allowing curtailment of renewable generators ``False``
``'optsolver'``             OPTALG optimization solver ``{'iqp','augl','ipopt'}`` ``'iqp'``
=========================== ===================================================== =========

The following example illustrates how to solve a DC OPF problem and extract the optimal generation cost::

  >>> method = gridopt.power_flow.new_method('DCOPF')

  >>> method.solve(net)

  >>> print method.results['status']
  solved

  >>> method.update_network(net)

  >>> # generation cost ($/hour)
  >>> print net.gen_P_cost
  7643.04

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
  gen  2 3.23e-02 1.02e+02

As the examples show, GRIDOPT and `PFNET`_ take care of all the details and allow one to extract solution information easily and intuitively from the network components.

.. _ac_pf:

ACPF
====

This method is represented by an object of type :class:`ACPF <gridopt.power_flow.ac_pf.ACPF>` and solves an AC power flow problem. For doing this, it can use the `Newton-Raphson`_ solver from `OPTALG`_ together with "switching" heuristics for modeling local controls. Alternatively, it can formulate the problem as an optimization problem with a convex objective function and `complementarity constraints`_ for modeling local controls, and solve it using the `Augmented Lagrangian`_, `INLP`_, or `IPOPT`_ solver from `OPTALG`_. For now, the parameters of this power flow method are the following:

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
``'optsolver'``   OPTALG optimization solver ``{'nr','inlp','augl','ipopt'}``  ``'augl'``
================= ============================================================ ===========

.. _ac_opf: 

ACOPF
=====

This method is represented by an object of type :class:`ACOPF <gridopt.power_flow.ac_opf.ACOPF>` and solves an AC optimal power flow problem. For doing this, it uses the `Augmented Lagrangian`_, `INLP`_, or `IPOPT`_ solver from `OPTALG`_. By default, it minimizes `active power generation cost`_ subject to voltage magnitude limits, generator power limits (`variable limits`_), and `AC power balance constraints`_. For now, the parameters of this optimal power flow method are the following:

==================== ============================================================ ===========
Name                 Description                                                  Default  
==================== ============================================================ ===========
``'weight_cost'``    Weight for active power generation cost                      ``1e0`` 
``'weight_mag_reg'`` Weight for soft voltage magnitude limits (or regularization) ``0.``
``'weight_ang_reg'`` Weight for voltage angle regularization                      ``0.``
``'weight_gen_reg'`` Weight for generator powers regularization                   ``0.``
``'thermal_limits'`` Flag for thermal limits                                      ``False``
``'vmin_thresh'``    Low-voltage threshold                                        ``1e-1``
``'optsolver'``      OPTALG optimization solver ``{'augl','inlp','ipopt'}``       ``'augl'``
==================== ============================================================ ===========

.. _OPTALG: http://optalg.readthedocs.io
.. _IQP: http://optalg.readthedocs.io/en/latest/opt_solver.html#iqp
.. _INLP: http://optalg.readthedocs.io/en/latest/opt_solver.html#inlp
.. _newton-raphson: http://optalg.readthedocs.io/en/latest/opt_solver.html#nr
.. _augmented lagrangian: http://optalg.readthedocs.io/en/latest/opt_solver.html#augl
.. _ipopt: http://optalg.readthedocs.io/en/latest/opt_solver.html#ipopt

.. _PFNET: http://pfnet-python.readthedocs.io
.. _network: http://pfnet-python.readthedocs.io/en/latest/networks.html
.. _contingencies: http://pfnet-python.readthedocs.io/en/latest/networks.html#contingencies
.. _optimization problem: http://pfnet-python.readthedocs.io/en/latest/problems.html#problems
.. _DC power balance constraints: http://pfnet-python.readthedocs.io/en/latest/problems.html#dc-power-balance
.. _AC power balance constraints: http://pfnet-python.readthedocs.io/en/latest/problems.html#ac-power-balance
.. _DC branch flow limits: http://pfnet-python.readthedocs.io/en/latest/problems.html#dc-branch-flow-limits
.. _variable limits: http://pfnet-python.readthedocs.io/en/latest/problems.html#variable-bounds
.. _active power generation cost: http://pfnet-python.readthedocs.io/en/latest/problems.html#active-power-generation-cost
.. _active power consumption utility: http://pfnet-python.readthedocs.io/en/latest/problems.html#active-power-consumption-utility
.. _network properties: http://pfnet-python.readthedocs.io/en/latest/networks.html#properties
.. _complementarity constraints: http://pfnet-python.readthedocs.io/en/latest/problems.html#voltage-set-point-regulation-by-generators

