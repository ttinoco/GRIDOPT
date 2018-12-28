.. include:: defs.hrst

.. _start:

***************
Getting Started
***************

This section describes how to get started with GRIDOPT. In particular, it covers installation and provides a quick example that shows how to use this package.

.. _start_installation:

Installation
============

In order to install GRIDOPT, the following tools are needed:

* Linux and macOS:

  * C compiler
  * |make|
  * |python| (2 or 3)
  * |pip|
  
* Windows:

  * |anaconda| (for Python 2.7)
  * |cmake| (choose "Add CMake to the system PATH for all users" during installation)
  * |7-zip| (update system path to include the 7z executable, typically in ``C:\Program Files\7-Zip``)
  * |mingwpy| (use ``pip install -i https://pypi.anaconda.org/carlkl/simple mingwpy``)

After getting these tools, the GRIDOPT Python module can be easily installed by executing the following commands on the terminal or Anaconda prompt::

  pip install numpy cython
  pip install optalg pfnet
  pip install gridopt

To install the module from source, the code can be obtained from `<https://github.com/ttinoco/GRIDOPT>`_, and then the following commands can be executed on the terminal or Anaconda prompt from the root directory of the package::

    pip install numpy cython
    pip install optalg pfnet
    python setup.py install

Running the unit tests can be done with::

    pip install nose
    nosetests -s -v

.. _start_example:

Example
=======

The following example shows how to solve the power flow problem associated with a power grid using GRIDOPT::

  >>> import pfnet
  >>> import gridopt

  >>> net = pfnet.PyParserMAT().parse('ieee14.m')

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
  5.16e-04 5.67e-03

In this example, it is assumed that the Python interpreter is started from a directory where the sample case |ieee14| is located.
