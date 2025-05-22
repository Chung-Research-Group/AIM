.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AIM Modules
===============================

**AIM** consists of 4 modules:

* **IsoFit**
* **HeatFit**
* **MixPred**
* **BreakLab**

.. image:: images/AIM_modules.png
   :width: 1000
   :alt: AIM_logo
   :align: center
   
-------
 
| **IsoFit**: Fit various isotherm models to single temperature isotherm data.
| **HeatFit**: Fit various isotherm models to multi-temperature isotherm data and predict isosteric heat of adsorption using Clausius-Clapeyron or Virial equation.
| **MixPred**: Predict mixture isotherm using Ideal Adsorbed Solution Theory (IAST) and Extended Dual-site Langmuir (EDSL) models.
| **BreakLab**: Predict multicomponent isothermal/nonisothermal breakthrough curves for fixed bed.



The isotherm fitting results from **IsoFit** and **HeatFit** can be directly loaded into **MixPred** and **BreakLab** modules,
providing seamless integration.

.. toctree::
   :maxdepth: 2
   :caption: Modules:

   IsoFit
   HeatFit
   MixPred
   BreakLab
