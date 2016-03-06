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

.. autoclass:: gridopt.power_flow.DCPF

.. autoclass:: gridopt.power_flow.DCOPF

.. autoclass:: gridopt.power_flow.NRPF

.. autoclass:: gridopt.power_flow.AugLPF

.. autoclass:: gridopt.power_flow.AugLOPF

.. autoclass:: gridopt.power_flow.method_error.PFmethodError
   :members:

.. _ref_pf_error:

Error Exceptions
----------------

.. autoclass:: gridopt.power_flow.method_error.PFmethodError

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_NoProblem

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadProblem

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadFlowLimits

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadVarLimits

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_BadParam

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_ParamNotBool

.. autoclass:: gridopt.power_flow.method_error.PFmethodError_SolverError
  
.. _ref_references:

References
==========

.. [TTR2015] T\. Tinoco De Rubira, *Numerical Optimization and Modeling Techniques for Power System Operations and Planning*. PhD thesis, Stanford University, March 2015.



