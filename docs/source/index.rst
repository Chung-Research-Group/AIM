AIM — Adsorption Integrated Modules
====================================

.. rst-class:: lead

**AIM** is an open-source MATLAB application suite that connects pure-component
isotherm fitting, isosteric heat estimation, mixture prediction, and fixed-bed
breakthrough simulation in one reproducible workflow.

.. grid:: 2 2 4 4
   :gutter: 3
   :class-container: aim-module-grid

   .. grid-item-card:: IsoFit
      :link: isofit
      :link-type: doc

      Fit single-temperature adsorption isotherms with 13 model families,
      bounded nonlinear regression, multistart search, and automatic model
      selection.

   .. grid-item-card:: HeatFit
      :link: heatfit
      :link-type: doc

      Fit multi-temperature data and estimate isosteric heat using the
      Clausius–Clapeyron or virial formulation.

   .. grid-item-card:: MixPred
      :link: mixpred
      :link-type: doc

      Predict multicomponent equilibrium with extended dual-site Langmuir
      (EDSL) or ideal adsorbed solution theory (IAST).

   .. grid-item-card:: BreakLab
      :link: breaklab
      :link-type: doc

      Simulate nonisothermal, non-isobaric fixed-bed breakthrough for as many
      as five components.

Start here
----------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Install AIM
      :link: installation
      :link-type: doc

      Use the Windows standalone installer or run the application from MATLAB
      R2024a or newer.

   .. grid-item-card:: Watch the tutorials
      :link: tutorials
      :link-type: doc

      Follow complete video walkthroughs for IsoFit, HeatFit, MixPred, and
      BreakLab.

   .. grid-item-card:: Prepare input files
      :link: input-output
      :link-type: doc

      Review supported formats, metadata tags, units, and the human-readable
      ``.bliso`` output contract.

   .. grid-item-card:: Reproduce published cases
      :link: examples
      :link-type: doc

      Use the datasets, parameter files, and BreakLab configurations under
      ``manuscript_data``.

Scientific scope
----------------

AIM retains the model definitions, governing equations, boundary conditions,
parameter units, fitting objectives, and numerical-method descriptions needed
to understand and reproduce its calculations. The module guides are technical
references, not only interface tours.

.. toctree::
   :hidden:
   :maxdepth: 3

   installation
   tutorials
   input-output
   modules
   isofit
   heatfit
   mixpred
   breaklab
   examples
   development
   citation

Project links
-------------

* `Source repository <https://github.com/Chung-Research-Group/AIM>`_
* `Issue tracker <https://github.com/Chung-Research-Group/AIM/issues>`_
* `Associated article <https://doi.org/10.1016/j.cpc.2025.109944>`_
