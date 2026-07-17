AIM modules
===========

IsoFit
------

Fits pure-component adsorption data using models including single- and
dual-site Langmuir, Langmuir–Freundlich, quadratic, Temkin, BET, Sips, Toth,
structural-transition, Dubinin–Astakhov, Klotz, and Do–Do formulations.

HeatFit
-------

Uses isotherms measured at multiple temperatures to estimate isosteric heat of
adsorption with Clausius–Clapeyron or virial methods.

MixPred
-------

Predicts multicomponent equilibrium adsorption with the extended dual-site
Langmuir model or ideal adsorbed solution theory (IAST).

BreakLab
--------

Simulates nonisothermal, non-isobaric fixed-bed breakthrough for as many as
five components. The model supports axial dispersion, linear-driving-force
mass transfer, and Ergun pressure drop.

Recommended workflow
--------------------

1. Fit pure-component data in **IsoFit**.
2. Estimate temperature dependence in **HeatFit**, when multi-temperature data
   are available.
3. Check mixture equilibrium predictions in **MixPred**.
4. Supply equilibrium and process parameters to **BreakLab**.
5. Export inputs and results with the case so the calculation can be
   reproduced.
