.. _start:

***************
Getting Started
***************

This section describes how to get started with GRIDOPT. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

* `Numpy <http://www.numpy.org>`_ (>=1.11.2)
* `Scipy <http://www.scipy.org>`_ (>=0.18.1)
* `OPTALG <http://optalg.readthedocs.io>`_ (>= 1.1.2)
* `PFNET`_ (>= 1.2.9)

.. _start_download:

Download
========

The latest version of GRIDOPT can be downloaded from `<https://github.com/ttinoco/GRIDOPT>`_.

.. _start_installation:

Installation
============

The GRIDOPT Python module can be installed using::

  > sudo python setup.py install

from the root directory of the package.

The installation can be tested using `nose <https://nose.readthedocs.org/en/latest/>`_ as follows::

  > nosetests -v -s

.. _start_example:

Example
=======

The next example shows how to solve the power flow problem associated with a power grid using GRIDOPT::

  >>> import pfnet
  >>> import gridopt

  >>> net = pfnet.ParserMAT().parse('ieee14.mat')

  >>> # max mismatches (MW,MVAr)
  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  3.54e-01 4.22e+00

  >>> method = gridopt.power_flow.new_method('NRPF')

  >>> method.set_parameters({'quiet': True})

  >>> method.solve(net)

  >>> results = method.get_results()

  >>> print results['status']
  solved

  >>> method.update_network(net)

  >>> # max mismatches (MW,MVAr)
  >>> print '%.2e %.2e' %(net.bus_P_mis,net.bus_Q_mis)
  5.16e-04 5.67e-03

In this example, it is assumed that the Python interpreter was started from the directory ``tests/resources/cases`` of the GRIDOPT package, where the sample case ``ieee14.mat`` is located.

.. _PFNET: http://pfnet-python.readthedocs.io/
