.. include:: defs.hrst

.. _start:

***************
Getting Started
***************

This section describes how to get started with GRIDOPT. In particular, it covers dependencies, installation, and provides a quick example showing how to use this package.

.. _start_dependencies:

Dependencies
============

GRIDOPT has the following dependencies:

* |Cython| (>=0.20.1)
* |Numpy| (>=1.11.2)
* |Scipy| (>=0.18.1)
* |OPTALG| (>= 1.1.5r1)
* |PFNET| (>= 1.3.2r1)

.. _start_installation:

Installation
============

In order to install the GRIDOPT, the following tools are needed:

* Linux and Mac OS X: a C compiler, |Make|, |Python| and |pip|.
* Windows : |Anaconda|, |CMake|, |7-Zip|, and |MinGW|.

After getting these tools, the GRIDOPT Python module can be easily installed by executing the following commands on the terminal or Anaconda prompt::

  pip install numpy cython
  pip install optalg pfnet
  pip install gridopt

To install the module from source, the code can be obtained from `<https://github.com/ttinoco/GRIDOPT>`_, and then the following commands can be executed on the terminal or Anaconda prompt from the root directory of the package::

    pip install numpy cython
    pip install optalg pfnet
    python setup.py install

Running the unit tests can be done with::

    nosetests -s -v

.. _start_example:

Example
=======

The following example shows how to solve the power flow problem associated with a power grid using GRIDOPT::

  >>> import pfnet
  >>> import gridopt

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> # max mismatches (MW,MVAr)
  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  3.54e-01 4.22e+00

  >>> method = gridopt.power_flow.new_method('ACPF')

  >>> method.set_parameters({'quiet': True})

  >>> method.solve(net)

  >>> results = method.get_results()

  >>> print results['solver status']
  solved

  >>> method.update_network(net)

  >>> # max mismatches (MW,MVAr)
  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.14e-04 5.70e-03

In this example, it is assumed that the Python interpreter was started from a directory where the sample case |ieee14| is located.
