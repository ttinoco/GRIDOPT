.. _script:

********************
Command-Line Utility
********************

This section provides details about how to use the ``gridopt`` command-line utility.

.. _script_syntax:

Syntax
======

The command-line utility ``gridopt`` can be used to load a network data file, and solve the corresponding PF or OPF problem using an available methods. One can also configure the parameters of the underlying optimization solver used as well as the properties of the optimization problem constructed. Saving solved cases to a file is not yet supported. The syntax and a description of each of the options is presented below.

.. program:: gridopt_clu

::

   usage: gridopt case 
                  method 
                  [--params <name1=value1> <name2=value2> ...] 
                  [--profile]
                  [--flatstart]

.. option:: case 

            Name of a power network data file.

.. option:: method

	    Name of method (``DCPF``, ``DCOPF``, ``ACPF``, ``ACOPF``).

.. option:: --params <name1=value1> <name2=value2> ...

	    List of parameter name-value pairs for configuring the method and the underlying optimization solver.

.. option:: --profile
	   
	    Profiles method execution using `cProfile <http://docs.python.org/2/library/profile.html#module-cProfile>`_.

.. option:: --flatstart

	    Enforces flat starting point (zero phase angles and unity voltage manigtudes).

.. _script_example:

Example
=======

The following example shows how to use the command-line utility to solve an AC power flow problem using the Newton-Raphson algorithm with a feasibility tolerance of ``1e-5`` per unit system MVA::

  gridopt ieee14.mat ACPF --params feastol=1e-5 solver=nr
