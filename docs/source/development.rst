Development and validation
==========================

Requirements
------------

* MATLAB R2024a or newer
* MATLAB Test and Optimization Toolbox for the CI test configuration
* MATLAB Compiler only when packaging a standalone application

Checks
------

From the repository root, run:

.. code-block:: matlab

   buildtool check
   buildtool test

Pull requests run repository hygiene, MATLAB Code Analyzer, and numerical
tests. The test matrix covers R2024a on Linux and the latest MATLAB release on
Windows. The required checks and artifact policy are documented in
`.github/CI.md <https://github.com/Chung-Research-Group/AIM/blob/main/.github/CI.md>`_.

Build these docs locally
------------------------

.. code-block:: console

   python -m pip install -r docs/requirements.txt
   sphinx-build -W --keep-going -b html docs/source docs/_build/html

Read the Docs uses the same pinned dependencies and treats warnings as build
failures. Its configuration is stored in ``.readthedocs.yaml`` at the
repository root.

Contributing a fix
------------------

Keep scientific changes small enough to review independently. Add a regression
test that fails before the fix, state the governing invariant or expected
result, and report any numerical or performance impact in the pull request.
