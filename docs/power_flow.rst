.. _power_flow:

**********
Power Flow
**********

The Power Flow (PF) problem consists of determining steady-state voltage magnitudes and angles at every bus of the network as well as any unknown generator powers. On the other hand, the Optimal Power Flow (OPF) problem consists of determining generator and other network control settings that result in the optimal operation of the network according to some measure, *e.g.*, generation cost. In GRIDOPT, methods for solving PF and OPF problems are represented by objects derived from a :class:`method <gridopt.power_flow.method.PFmethod>` base class.

To solve a PF or OPF problem, one first needs to create an instance of a specific method subclass. This is done using the function :func:`new_method() <gridopt.power_flow.new_method>`, which takes as argument the name of an available power flow method (:ref:`dc_pf`, :ref:`dc_opf`, :ref:`nr_pf`, :ref:`augl_pf`, and :ref:`augl_opf`). The following code sample creates an instance of the :ref:`nr_pf` method::

  >>> import gridopt

  >>> method = gridopt.power_flow.new_method('NRPF')

Once a method has been instantiated, its parameters can be set using the function :func:`set_parameters() <gridopt.power_flow.method.PFmethod.set_parameters>`. This function takes as argument a dictionary with parameter name and value pairs. Valid parameters include parameters of the method, which are described in the section below, and parameters from the numerical solver used by the method. The numerical solvers used by the methods of GRIDOPT belong to the Python package `OPTALG <https://github.com/ttinoco/OPTALG>`_. The following code sample sets a few parameters of the method created above::

  >>> method.set_parameters({'quiet': True, 'feastol': 1e-4})

After configuring parameters, a method can be used to solve a problem using the function :func:`solve() <gridopt.power_flow.method.PFmethod.solve>`. This function takes as argument a `PFNET <http://ttinoco.github.io/PFNET/python/>`_ Network object, as follows::

  >>> import pfnet

  >>> net = pfnet.Network()
  >>> net.load('ieee14.mat')

  >>> method.solve(net)

Information about the execution of the method can be obtained from the :data:`results <gridopt.power_flow.method.PFmethod.results>` attribute of the :class:`method <gridopt.power_flow.PFmethod.results>` object. This dictionary of results includes information such as ``status``, *e.g.*, ``solved`` or ``error``, any error message (``error_msg``), solver ``iterations``, and `PFNET <http://ttinoco.github.io/PFNET/python/>`_ optimization ``problem`` constructed, and the properties of the `PFNET <http://ttinoco.github.io/PFNET/python/>`_ Network associated with the point found by the method (``net_properties``). The following code sample show how to extract and see some results::

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

If desired, one can update the `PFNET <http://ttinoco.github.io/PFNET/python/>`_ Network object with the solution found by the method. This can be done with the function :func:`update_network() <gridopt.power_flow.method.PFmethod.update_network>`. This routine not only updates the network quantities treated as variables by the method, but also information about sensitivities of the objective function (if any) with respect to perturbation of the constraints. The following code sample update the power network with the results obtained by the method and shows the resulting maximum active and reactive bus power mismatches in units of MW and MVAr::

  >>> method.update_network(net)

  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.16e-04 5.67e-03
    
.. _dc_pf: 

DCPF
====

.. _dc_opf: 

DCOPF
=====

.. _nr_pf: 

NRPF
====

.. _augl_pf: 

AugLPF
======

This method is similar to the algorithm described in Chapter 3 of [TTR2015]_, but without the restriction of moving in the null-space of the linear equality constraints.

.. _augl_opf: 

AugLOPF
=======


