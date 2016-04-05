.. _power_flow:

**********
Power Flow
**********

The Power Flow (PF) problem consists of determining steady-state voltage magnitudes and angles at every bus of the network as well as any unknown generator powers. On the other hand, the Optimal Power Flow (OPF) problem consists of determining generator and other network control settings that result in the optimal operation of the network according to some measure, *e.g.*, generation cost. In GRIDOPT, methods for solving PF and OPF problems are represented by objects derived from a :class:`method <gridopt.power_flow.method.PFmethod>` base class.

To solve a PF or OPF problem, one first needs to create an instance of a specific method subclass. This is done using the function :func:`new_method() <gridopt.power_flow.new_method>`, which takes as argument the name of an available power flow method. Available power flow methods are the following: 

* :ref:`dc_pf`
* :ref:`dc_opf`
* :ref:`dc_opf_mp`
* :ref:`nr_pf`
* :ref:`augl_pf`
* :ref:`augl_opf`.

The following code sample creates an instance of the :ref:`nr_pf` method::

  >>> import gridopt

  >>> method = gridopt.power_flow.new_method('NRPF')

Once a method has been instantiated, its parameters can be set using the function :func:`set_parameters() <gridopt.power_flow.method.PFmethod.set_parameters>`. This function takes as argument a dictionary with parameter name and value pairs. Valid parameters include parameters of the method, which are described in the sections below, and parameters from the numerical solver used by the method. The numerical solvers used by the methods of GRIDOPT belong to the Python package `OPTALG`_. The following code sample sets a few parameters of the method created above::

  >>> method.set_parameters({'quiet': True, 'feastol': 1e-4})

After configuring parameters, a method can be used to solve a problem using the function :func:`solve() <gridopt.power_flow.method.PFmethod.solve>`. Typically, this function takes as argument a `PFNET Network`_ object, as follows::

  >>> import pfnet

  >>> net = pfnet.Network()
  >>> net.load('ieee14.mat')

  >>> method.solve(net)

Information about the execution of the method can be obtained from the :data:`results <gridopt.power_flow.method.PFmethod.results>` attribute of the :class:`method <gridopt.power_flow.method.PFmethod>` object. This dictionary of results includes information such as ``'status'``, *e.g.*, ``'solved'`` or ``'error'``, any error message (``'error_msg'``), solver ``'iterations'``, the `PFNET Optimization Problem`_ (``'problem'``) constructed, and the properties of the `PFNET Network`_ at the point found by the method (``'net_properties'``). The following code sample show how to extract some results::

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
    type: FIX
    type: PAR_GEN_Q
    type: PAR_GEN_P
    type: PF

  >>> print results['net_properties']['bus_v_max']
  1.09

If desired, one can update the `PFNET Network`_ object with the solution found by the method. This can be done with the function :func:`update_network() <gridopt.power_flow.method.PFmethod.update_network>`. This routine not only updates the network quantities treated as variables by the method, but also information about the sensitivity of the optimal objective function value with respect to perturbations of the constraints. The following code sample updates the power network with the results obtained by the method and shows the resulting maximum active and reactive bus power mismatches in units of MW and MVAr::

  >>> method.update_network(net)

  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.16e-04 5.67e-03
    
.. _dc_pf: 

DCPF
====

This method solves a DC power flow problem, which is just a linear system of equations. For doing this, it uses a `PFNET DC power balance constraint`_ and one of the ``linear solvers`` from `OPTALG`_.

.. _dc_opf: 

DCOPF
=====

This method solves a DC optimal power flow problem, which is just a quadratic program including active power generation cost, active power consumption utility, power balance, generator limits, load limits, and branch thermal limits. For doing this, it uses the ``OptSolverIQP`` interior point solver from `OPTALG`_. 

The following example illustrates how to solve a DCOPF problem and extract the optimal generation cost::

  >>> method = gridopt.power_flow.new_method('DCOPF')

  >>> method.solve(net)

  >>> print method.results['status']
  solved

  >>> method.update_network(net)

  >>> # generation cost ($/hour)
  >>> print net.gen_P_cost
  4810.98

The sensitivity of the optimal generation cost with respect to the power balance equations can be easily extracted from the network buses::

  >>> bus = net.get_bus(4)
  >>> print "bus %2d %.2e" %(bus.index,bus.sens_P_balance)
  bus 4 2.13e+01
  
Similarly, the sensitivity with respect to branch flow limits can be easily extracted from the network branches::

  >>> branch = net.get_branch(6)
  >>> print "branch %2d %.2e %.2e" %(branch.index,
  ...                                branch.sens_P_u_bound,
  ...                                branch.sens_P_l_bound)
  branch 6 2.01e-09 1.25e-09

Lastly, the sensitivity with respect to generator active power limits can be easily extracted from the network generators::

  >>> gen = net.get_gen(2)
  >>> print "gen %2d %.2e %.2e" %(gen.index,
  ...                             gen.sens_P_u_bound,
  ...                             gen.sens_P_l_bound)
  gen  2 2.01e-06 2.85e+01

As the examples show, GRIDOPT and PFNET take care of all the details and allow one to extract solution information easily and intuitively from the network components.

.. _dc_opf_mp: 

DCOPF_MP
========

This method solves a multi-period DC optimal power flow problem.

.. _nr_pf: 

NRPF
====

This method solves an AC power flow problem, which is a nonlinear system of equations. For doing this, it uses the ``OptSolverNR`` Newton-Raphson solver from `OPTALG`_. For now, its parameters are a ``'quiet'`` flag and a low-voltage threshold ``'vmin_thresh'``.

.. _augl_pf: 

AugLPF
======

This method solves an AC power flow problem but formulated as an optimization problem with a strongly-convex objective function. For doing this, it uses the ``OptSolverAugL`` Augmented Lagrangian solver from `OPTALG`_. The ``OptSolverAugL`` solver is similar to the one described in Chapter 3 of [TTR2015]_, but without the restriction of moving in the null-space of the linear equality constraints. For now, the parameters of this power flow method are the following:

================= ================================================ ===========
Name              Description                                      Default  
================= ================================================ ===========
``'weight_vmag'`` Weight for bus voltage magnitude regularization  ``1e0``
``'weight_vang'`` Weight for bus voltage angle regularization      ``1e-3``
``'weight_pq'``   Weight for generator power regularization        ``1e-3``
``'weight_t'``    Weight for transformer tap ratio regularization  ``1e1``
``'weight_b'``    Weight for shunt susceptance regularization      ``1e-4``
``'vmin_thresh'`` Low-voltage threshold                            ``1e-1``
================= ================================================ ===========

.. _augl_opf: 

AugLOPF
=======

This method solves an AC optimal power flow problem. For doing this, it uses the ``OptSolverAugL`` Augmented Lagrangian solver from `OPTALG`_. For now, the parameters of this optimal power flow method are the following:

================== ================================================ ===========
Name               Description                                      Default  
================== ================================================ ===========
``'weight_cost'``  Weight for active power generation cost          ``1e-2`` 
``'weight_limit'`` Weight for soft constraint violations            ``1e-2``
``'weight_reg'``   Weight for regularization                        ``1e-5``
``'vmin_thresh'``  Low-voltage threshold                            ``1e-1``
================== ================================================ ===========

.. _OPTALG: http://ttinoco.github.io/OPTALG/
.. _PFNET Network: http://ttinoco.github.io/PFNET/python/reference.html#network
.. _PFNET Optimization Problem: http://ttinoco.github.io/PFNET/python/reference.html#optimization-problem
.. _PFNET DC power balance constraint: http://ttinoco.github.io/PFNET/python/problems.html#dc-power-balance
