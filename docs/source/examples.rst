Example and publication data
============================

The ``manuscript_data`` directory contains inputs used for case studies in the
associated AIM article. It includes CSV isotherm datasets, saved ``.bliso``
parameter files, mixture-prediction data, and BreakLab configuration files.

To reproduce a case:

1. Clone or download the complete repository so companion files stay together.
2. Select the module named by the case-study directory.
3. Import the provided data or configuration through that module's interface.
4. Record the AIM commit, MATLAB or Runtime version, operating system, and any
   changed parameters with exported results.

Browse the
`manuscript data <https://github.com/Chung-Research-Group/AIM/tree/main/manuscript_data>`_
or read the
`data README <https://github.com/Chung-Research-Group/AIM/blob/main/manuscript_data/README.md>`_.

Reproducibility checklist
-------------------------

* Keep raw inputs separate from exported or transformed files.
* Preserve units and component order exactly as entered.
* Save fitted isotherm parameters used by MixPred or BreakLab.
* Report solver warnings and convergence failures; do not silently discard
  them.
