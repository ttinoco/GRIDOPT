Targets
-------

Unreleased
----------
* Reduced width of ACOPF output per iteration (done).
* Updated ACOPF regularization param names to be more consistent with those of ACPF (done).
* Hid internal param data structure from user in power flow methods, and allowed saving/restoring params including optsolver parameters (done).
* Changed method results and mechanism for storing solution and updating network (done).
* Moved location of command-line unitlity (done).
* Fixed bug in NR ACPF related to tap_changers vs tap_changers_v (done).
* Changed travis to pfnet from testpypi (done).
* Added "solver name" to method results (done).
* Changed "optsolver" to "solver" and added "solver" to DCPF (done).
* Add unittests for command-line utility (todo).

Version 1.3.3
-------------
* Unified methods (ACPF, ACOPF, DCPF, DCOPF) with param "optsolver" for choosing optimization solver.
* Added documentation for "gridopt" command-line utility.
* Updated and re-enabled DCOPF unittests.
* Added new OPTALG solver "inlp" to ACPF and ACOPF and updated documentation.
* limit_gens flag for non-NR ACPF with unittest.

Version 1.3.2
-------------
* Compatibility with PFNET v1.3.0.
* Compatibility with OPTALG v1.1.3.
* Thermal limits option for AugL and Ipopt OPFs.

Version 1.3.1
-------------
* Compatibility with PFNET 1.2.9.
* Tool for transforming PFNET to OPTALG problems and recovering primal and dual variables.
* Standard OPF problem formulation.
* OPF method based on AugL uses linear variable bounds and handles them with barrier.
* OPF method based on IPOPT solver.
* Flat start option in gridopt script.
* Travis continuous integration.
* Readthedocs.

Version 1.3
-----------
* Compatibility with multiperiod pfnet (dc_opf, nr_pf, augL_pf and augL_opf)

Version 1.2
-----------
* Multi-stage stochastic DCOPF.
* Variable generator power curtailment in DCOPF.
* Python 2 and 3 compatibility.

Version 1.1
-----------
* Multi-period DCOPF.
* Corrective DCOPF.
* Python 3 and Jupyter compatibility.
* Updated Dockerfile.
* Elastic loads in DCOPF.
* Thermal limits flag in DCOPF.

Version 1.0
-----------
* Initial version.
