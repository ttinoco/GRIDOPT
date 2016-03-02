.. _pf_methods:

******************
Power Flow Methods
******************

The power flow problem consists of determining steady-state voltage magnitudes and angles at every bus of the network as well as any unknown generator powers. In GRIDOPT, power flow methods are collectively represented by objects derived from a :class:`method <gridopt.power_flow.method>` base class.

To solve a power flow problem, one needs to first create an instance of a specific method subclass. This is done using the function :func:`new_method <gridopt.power_flow.new_method>`, which takes as argument the name of an available power flow method (``DCPF``, ``NRPF``, ``AugLPF``). These methods are described in the sections below. The following example creates an instance of the :ref:`pf_nr`::

  

.. _pf_dc: 

DC Power Flow Method
====================

.. _pf_nr:

Newton-Raphson Power Flow Method
================================

.. _pf_augl:

Augmented Lagrangian Power Flow Method
======================================

This method is similar to the algorithm described in Chapter 3 of [TTR2015]_, but without the restriction of moving in the null-space of the linear equality constraints.
