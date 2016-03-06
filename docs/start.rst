.. _start:

***************
Getting Started
***************

This section describes how to get started with GRIDOPT. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

* `Numpy <http://www.numpy.org>`_ (>=1.8.2)
* `Scipy <http://www.scipy.org>`_ (>=0.13.3)
* `OPTALG <https://github.com/ttinoco/OPTALG>`_
* `PFNET <http://ttinoco.github.io/PFNET/python/>`_ (>= 1.1)

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

  > nosetests -v

.. _start_docker:

Docker
======

If GRIDOPT was obtained as a `Docker <https://www.docker.com/>`_ image, say ``gridopt.tar``, then one needs to first install `Docker Engine <https://docs.docker.com/engine/installation/>`_. Then one can load the image using the command::

  > docker load -i gridopt.tar

and enter the application environment with::

  > docker run -i -t --entrypoint=/bin/bash gridopt

In the application environment, GRIDOPT and all its dependencies, *e.g.*, `PFNET <http://ttinoco.github.io/PFNET/python/>`_, are already installed and ready to go. There, one can navidate to the directory ``/gridopt/tests/resources`` and use the test cases available there to do the `PFNET <http://ttinoco.github.io/PFNET/python/>`_ and GRIDOPT tutorials. 

.. _start_example:

Example
=======

The next example shows how to solve the power flow problem associated with a power grid using GRIDOPT::

  >>> import pfnet
  >>> import gridopt

  >>> net = pfnet.Network()
  >>> net.load('ieee14.mat')

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
