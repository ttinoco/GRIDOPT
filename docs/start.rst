.. _start:

***************
Getting Started
***************

This section describes how to get started with GRIDOPT. In particular, it covers required packages, installation, and provides a quick example showing how to use this package.

.. _start_requirements:

* `Numpy <http://www.numpy.org>`_ (>=1.8.2)
* `Scipy <http://www.scipy.org>`_ (>=0.13.3)
* `OPTALG <http://ttinoco.github.io/OPTALG>`_ (== 1.0)
* `PFNET`_ (>= 1.3)

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

If GRIDOPT was obtained as a `Docker <https://www.docker.com/>`_ image, say `gridopt.tar <https://drive.google.com/file/d/0B8iMq4kjB3lIZFVSbHhRSGs5VUE/view?usp=sharing>`_, then one needs to first install `Docker Engine <https://docs.docker.com/engine/installation/>`_. Then one can load the image using the command::

  > docker load -i gridopt.tar

and enter the application environment with::

  > docker run -i -t --entrypoint=/bin/bash gridopt

In the application environment, GRIDOPT and all its dependencies, *e.g.*, `PFNET`_, are already installed and ready to go. There, one can navidate to the directory ``/gridopt/tests/resources`` and use the test cases available there to do the `PFNET`_ and GRIDOPT tutorials with Python. 

.. _start_docker_lin:

Graphics in Linux
-----------------

To display graphics within the application environment, the following command can be run for entering the application environment::

  > docker run -i -t --entrypoint=/bin/bash -e DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix:ro gridopt

Then, on the host machine the command ``xhost +`` can be used to enable access to your host machine's display and then ``xhost -`` to disable it after usage. Inside the application environment, the command ``xeyes`` can be used to check whether graphics are working.

.. _start_docker_win:

Graphics in Windows
-------------------

Displaying graphics on Windows involves a few more steps. First, `Xming <https://sourceforge.net/projects/xming/>`_, an X server for Windows, must be downloaded and installed. Then, the installed application ``XLaunch`` should be executed with the options ``Multiple windows``, ``Display number`` 0, ``Start no client``, ``Clipboard`` and ``No Access Control``. Once this is done, the application environment can be entered using::

  > docker run -i -t --entrypoint=/bin/bash -e DISPLAY=ip_address_taken_by_Xming:0 gridopt

where the IP address taken by Xming can be found in the file ``C:\Users\username\AppData\Local\Temp\Xming.0`` next to ``XdmcpRegisterConnection``. Again, graphics within the application environment can be tested using the command ``xeyes``.

.. _start_docker_mac:

Graphics in Mac
---------------

Coming soon.

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

.. _PFNET: http://ttinoco.github.io/PFNET/python/
