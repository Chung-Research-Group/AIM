.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MixPred
===============================

**MixPred** is the module for mixture isotherm prediction for the given pressure and composition.
**MixPred** GUI is shown below:

.. image:: images/MixPred.png
   :width: 1000
   :alt: MixPred
   :align: center
   
-------

The isotherm fitting results from **IsoFit** and **HeatFit** can be directly loaded into **MixPred** module,
and the user can specify the desired pressure range and mixture composition to calculate the mixture loadings. 

Watch how to use **MixPred** :ref:`here <MixPred-label>`

**MixPred** has two models for
mixture adsorption loading prediction which are as follows:

* Extended Dual-site Langmuir Model (EDSL)
* Ideal Adsorbed Solution Theory (IAST)


Extended Dual-site Langmuir Model (EDSL)
---------------------------------------------

The extended dual-site Langmuir (EDSL) model is a generalization of the Langmuir model for multi-component mixture adsorption. The model is expressed as,

.. math::
  q_{mix,i}^{*} = \frac{q_{sat,1,i}b_{1,i}P_{i}}{1+\sum_{1}^{N}b_{1,j}P_{j}}
            + \frac{q_{sat,2,i}b_{2,i}P_{i}}{1+\sum_{1}^{N}b_{2,j}P_{j}}

where :math:`N` is the number of components involved; :math:`q_{sat,1,i}`, :math:`q_{sat,2,i}`, :math:`b_{1,i}`, :math:`b_{2,i}`, and :math:`P_{i}` are
the single or pure component Langmuir isotherm parameters and partial pressure of component :math:`i`, respectively. The pure component Langmuir isotherm parameters for 
individual components can be obtained by fitting Langmuir isotherm model to the pure component isotherm data for the given component using either **IsoFit** or **HeatFit**. 
Note that the EDSL model is thermodynamically consistent only when the saturation capacities for each component are equal.

.. math:: 
   q_{sat, 1, 1}=q_{sat, 1, 2}=q_{sat, 1, 3}= ... q_{sat, 1, N} \\
   q_{sat, 2, 1}=q_{sat, 2, 2}=q_{sat, 2, 3}= ... q_{sat, 2, N}

.. note::
   The EDSL model in **MixPred** is applicable only when the isotherm of each individual component is described using the Langmuir isotherm model.

Ideal Adsorbed Solution Theory (IAST)
---------------------------------------------

IAST is a thermodynamic framework to calculate the mixture isotherms using pure component isotherm. The IAST is based on three fundamental assumptions:

* The surface area of the adsorbent is equally accessible to all adsorbates.
* The adsorbed phase is an ideal solution.
* The adsorbent is homogeneous.

The solution of the IAST involves solving non-linear equations consisting of a rediced grand potential. For the component :math:`i`, the reduced grand potential :math:`\psi_{i}^{*}` is given as,

.. math::
   \psi_{i}^{*} = \int_{0}^{P_{i}^{*}} \frac{q_{i}^{*}(P)}{P} \,dP

where :math:`P_{i}^{*}` is the fictitious pressure and :math:`q_{i}^{*}` is the equilibrium loading of component :math:`i` as obtained from pure component isotherm model. The fictitious pressure :math:`P_{i}^{*}` is the pressure
for component :math:`i` at which it exerts the same reduced grand potential as the other components in the mixture. :math:`P_{i}^{*}` is related to partial pressure :math:`P_{i}` of component :math:`i` given as,

.. math::
   P_{i} = x_{i}\,P_{i}^{*} \qquad \mathrm{for}\, i=1,2,3,\ldots, N.

where :math:`N` is the number of adsorbing components in the mixture and :math:`x_{i}` is the adsorbed mole fraction of component :math:`i`. Note that,

.. math::
   \sum_{i}\,x_{i} = 1

IAST states that adsorption equilibrium is achieved when the reduced grand potential of each component is same,

.. math::
   \psi_{1}^{*} = \psi_{2}^{*} = \psi_{3}^{*} = ... \psi_{N}^{*} 

In **MixPred** the IAST equations above are solved for :math:`2N` unknows, :math:`P_{i}^{*}` and :math:`x_{i}`. The mixture adsorption loading for each component are then calculated using

.. math::
   q_{tot} = \frac{1}{\sum_{i}^{N}\,\left(\frac{x_{i}}{q_{i}^{*}(P_{i}^{*})}\right)} \\
   q_{mix, i}^{*} = x_{i}\,q_{tot}

.. note::
   The IAST equations are highly non-linear and often extremely sensitive to the initial guess used. The solution can often be computationally expensive and time-consuming.