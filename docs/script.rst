.. include:: defs.hrst

.. _script:

********************
Command-Line Utility
********************

This section provides details about how to use the ``gridopt`` command-line utility.

.. _script_syntax:

Syntax
======

The command-line utility ``gridopt`` can be used to load a network data file and solve the corresponding PF or OPF problem using an available method. One can also configure the parameters of the underlying optimization solver as well as the properties of the optimization problem constructed. The syntax and a description of each of the options is presented below.

.. argparse::
   :module: gridopt.scripts.gridopt
   :func: create_parser
   :prog: gridopt

.. _script_example:

Example
=======

The following example shows how to use the command-line utility to solve an AC power flow problem using the Newton-Raphson algorithm with a feasibility tolerance of ``1e-5`` per unit system MVA::

  gridopt ieee14.m ACPF --parameters feastol=1e-5 solver=nr

In this example, it is assumed that the command is executed from a directory where the sample case |ieee14| is located.
