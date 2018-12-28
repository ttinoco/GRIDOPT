.. include:: defs.hrst

.. _reference:

*************
API Reference
*************

.. _ref_pf:

Power Flow Method
=================

.. autofunction:: gridopt.power_flow.new_method

.. autoclass:: gridopt.power_flow.method.PFmethod
   :members:

.. autoclass:: gridopt.power_flow.dc_pf.DCPF

.. autoclass:: gridopt.power_flow.dc_opf.DCOPF

.. autoclass:: gridopt.power_flow.ac_pf.ACPF

.. autoclass:: gridopt.power_flow.ac_opf.ACOPF

.. _ref_pf_error:

Error Exceptions
----------------

.. autoclass:: gridopt.power_flow.method_error.PFmethodError
   :members:

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_NoProblem

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadProblem

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadFlowLimits

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadVarLimits

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadParams

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadOptSolver			      

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_ParamNotBool
			      
.. autoclass:: gridopt.power_flow.method_error.PFmethodError_SolverError
