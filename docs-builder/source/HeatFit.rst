.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HeatFit
===============================

**HeatFit** is the module for isotherm model fitting to multi-temperature isotherm data and
isosteric heat of adsoprtion :math:`\Delta H_{ads}` prediction using Clausius-Clapeyron or Virial equations. **HeatFit** GUI is shown below:

.. image:: images/HeatFit.png
   :width: 1000
   :alt: HeatFit
   :align: center
   
-------

Watch how to use **HeatFit** :ref:`here <HeatFit-label>`

In **HeatFit**, isotherm fitting and prediction of the isosteric heat of adsorption can be carried out using either the Clausius-Clapeyron or Virial equations.
Fitting method varies based on the chosen model.

Clausius-Clapeyron Equation
---------------------------------------------

The Clausius-Clapeyron is a thermodynamic relationship which expresses the relationship between pressure and temperature for phase change processes.
For adsorption processes, the Clausius-Clapeyron equation is given as,

.. math::
  \left. \frac{\partial(\ln P)}{\partial T} \right|_q = \frac{-\Delta H_{ads}}{RT^{2}}

where :math:`R` is the general gas constant. For Clausius-Clapeyron implementation in **HeatFit** we assume that the :math:`\Delta H_{ads}` is independent of temperature and
pressure. The assumption yields the following form:

.. math::
  \ln\left(\frac{P_{ref}}{P}\right) = \frac{-\Delta H_{ads}}{R}\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)

where :math:`P_{ref}`, :math:`T_{ref}`, :math:`P`, and :math:`T` are the adsorption pressure and temperatures for the reference and given states, respectively.
Note that both reference and given state correspond to the same adosrption loading :math:`q`. Let :math:`f` be the function representing isotherm model then, 

.. math::
  q = f(T, P)=f(T_{ref}, P_{ref})

Assuming :math:`f_{ref}` be the isotherm model fitted at reference temperature, we have

.. math::
  q = f(T, P)=f_{ref}(P_{ref})

The value of :math:`P_{ref}` can be subsitutted using the Clausius-Clapeyron quation shown above such that,

.. math::
  q = f(T, P)=f_{ref}\left(P\times\exp{\left(\frac{-\Delta H_{ads}}{R}\left(\frac{1}{T}-\frac{1}{T_{ref}}\right)\right)}\right)

The utility of the above equation is that it allows us to predict the loading at different temperatures and pressures using the isotherm model fitted at reference condition.
Additionally, if we have pressure-loading data at different temperatures and an isotherm model fitted at the reference conditions we can also fit for :math:`\Delta H_{ads}`. 


**HeatFit** uses the methodology outlined by Ga *et al.* [#Ga]_ for isotherm fitting at reference condition and :math:`\Delta H_{ads}` prediction. The method involves following step :

* **Step 1**: Isotherm fitting at reference temperature
* **Step 2**: Fitting for :math:`\Theta` parameter
* **Step 3**: Fitting for isosteric heat of adsorption :math:`\Delta H_{ads}`

In all three steps, regression is performed using MATLAB's built-in lsqnonlin solver.
The user can control the regression process using different sets of values for initial guesses, lower and upper bounds. 

Step 1: Isotherm Fitting at Reference Temperature
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like **IsoFit**, the isotherm fitting at reference temperature is performed for the chosen isotherm model using non-linear regression.
Currently, **HeatFit** supports following isotherm models:

.. list-table:: HeatFit Isotherm Models
   :header-rows: 1
   :widths: auto

   * - Isotherm Models
     - Isotherm Expression
     - Parameters
   * - Langmuir
     - :math:`q = \frac{q_{sat} b P}{1 + b P}`
     - :math:`q_{sat},\ b`
   * - Dual-site Langmuir
     - :math:`q = \frac{q_{sat,1} b_{1} P}{1 + b_{1} P} + \frac{q_{sat,2} b_{2} P}{1 + b_{2} P}`
     - :math:`q_{sat,1},\ b_{1},\ q_{sat,2},\ b_{2}`
   * - Langmuir-Freundlich
     - :math:`q = \frac{q_{sat} b P^{n}}{1 + b P^{n}}`
     - :math:`q_{sat},\ b,\ n`
   * - Dual-site Langmuir-Freundlich
     - :math:`q = \frac{q_{sat, 1} b P^{n_{1}}}{1 + b_{1} P^{n_{1}}} + \frac{q_{sat,2} b_{2} P^{n_{2}}}{1 + b_{2} P^{n_{2}}}`
     - :math:`q_{sat,1},\ b_{1},\ n_{1},\ q_{sat,2},\ b_{2},\ n_{2}`
   * - Quadratic
     - :math:`q = q_{sat}\left(\frac{b P + c P^{2}}{1 + b P + c P^{2}}\right)`
     - :math:`q_{sat},\ b,\ c`
   * - Temkin
     - :math:`q = q_{sat}\left(\frac{b P}{1 + b P}\right) + q_{sat}\theta\left(\frac{b P}{1 + b P}\right)\left(\frac{b P}{1 + b P}-1\right)`
     - :math:`q_{sat},\ b,\ \theta` 
   * - BET
     - :math:`q = \frac{q_{sat} b P}{(1 - c P)(1 - c P + b P)}`
     - :math:`q_{sat},\ b,\ c`
   * - Sips
     - :math:`q = \frac{q_{sat} (b P)^{1/n}}{1 + (b P)^{1/n}}`
     - :math:`q_{sat},\ b,\ n`
   * - Toth
     - :math:`q = \frac{q_{sat} b P}{\left(1 + (b P)^n\right)^{1/n}}`
     - :math:`q_{sat},\ b,\ n`


The objective of the non-linear regression is to minimize the :math:`SSE` function given as,

.. math::
  SSE = \min_{a_{k}} \sum_{i}^N \left(q_{i, exp}^{*} - f_{ref}(P_{i}; \{a_{k}\}_{1}^{M})\right)^{2}

where :math:`N` and :math:`q_{i,exp}^{*}` represents the total number of data points and the experimental gas uptake for the given data point :math:`i`, respectively. 
:math:`f_{ref}(P_{i}; \{a_{k}\}_{1}^{M})` is the isotherm model corresponding to reference conditions.
where, :math:`P_{i}` is the pressure value for the given data point :math:`i`, :math:`a_{k}` is the set of parameters
and :math:`M` is the total number of parameters for the given isotherm model.

In **HeatFit** the user can control the regression process by specifying custom initial guesses, as well as lower and upper bounds.
Additionally, **HeatFit** offers a multistart option, which generates 1000 random initial guess within the specified bounds.
The fitting process is then performed sequentially for each initial guess and the best fitting result is selected. The multistart approach is useful for fitting problems with multiple parameter solutions. In such cases, fitting using multistart option can identify the global minimum corresponding to the best parameter estimates. The multistart option is available for all isotherm models except **Auto** mode. In **Auto** mode, IsoFit performs isotherm fitting using all the available isotherm models and then selects the best model. Using multistart option for **Auto** mode can be computationally expensive leading to excessive running times.

Root mean square error :math:`(RMSE)` is used in the program to evaluate the goodness of fit.

.. math::
  RMSE = \sqrt{\frac{SSE}{N-M}}

If the user chooses the **Auto** mode, **HeatFit** reports the best isotherm model with the lowest RMSE values. 
If two or more models have same value of SSE, then **HeatFit** will choose the model with a smaller number of parameters because of lower value of RMSE for :math:`\Delta H_{ads}` predcition.

**HeatFit** also reports coefficient of determination, :math:`r^{2}`, value defined as:

.. math::
  r^{2} = 1 - \frac{SSE}{\sum (q_{i,exp} - \overline{q_{i,exp}})}

where :math:`\overline{q_{i,exp}}` is the mean value of experimental gas uptakes.

Step 2: Fitting for :math:`\Theta` Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Second step involves using the isotherm model and the parameters obtained from fitting at reference temperature to get a set of parameters :math:`\Theta`,
which is defined as,

.. math::
  \Theta = \{\theta_j\}_{j=1}^Z

where :math:`Z` is the total number of temperature values used in the fitting process.
Every :math:`θ_{j}` value in the :math:`\Theta` vector corresponds to temperature value :math:`T_{j}`.
Then the program uses non-linear regression to fit for each :math:`θ_{j}` value with the objective function defined as follows,

.. math::
  SSE = \min_{\theta_{j}} \left(q_{j, exp}^{*} - g(P; \theta_{j})\right)^{2}

where :math:`q_{j, exp}^{*}` is the experimental gas uptake at temperature :math:`T_{j}` while
:math:`g(P; \theta_{j})` can be expressed as,

.. math::
  g(P; \theta_{j}) = f_{ref} \left(P\times\theta_{j};\ \{a_{k,ref}\}_{k=1}^{M}\right)

where :math:`f_{ref}` and :math:`a_{k,ref}` represent the isotherm function and isotherm function
parameters obtained in step 1 at reference conditions, respectively.

Step 3: Fitting for :math:`\Delta H_{ads}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, the correlation between :math:`θ_{j}` and :math:`T_{j}` is used to fit for :math:`\Delta H_{ads}`. 
The correlation is given by the Clausius-Clapeyron equation as:

.. math::
  \theta_{j, pred} = \exp\left(\frac{-\Delta H_{ads}}{R}\left(\frac{1}{T_{j}}-\frac{1}{T_{ref}}\right)\right)

The program uses non-linear regression to fit for :math:`\Delta H_{ads}` value 
while minimizing the objective function defined as:

.. math::
  SSE = \min_{\Delta H_{ads}} \sum_{j=1}^{Z}\left(\theta_{j} - \theta_{j, pred}\right)^{2}

**HeatFit** also reports RMSE value for :math:`\Delta H_{ads}` fitting, which is given as,

.. math::
  RMSE = \sqrt{\frac{SSE}{Z}}

---------------------------------------------

Virial Equation
---------------------------------------------

The Virial equation is given as,

.. math::
  \ln P = \ln q^{*} + \frac{1}{T}\sum_{i=0}^{m-1}a_{i} (q^{*})^{i} + \sum_{j=0}^{n-1}b_{j} (q^{*})^{j}

where :math:`a_{i}` and :math:`b_{j}` are the Virial fitting parameters also known as Virial coefficients. :math:`m` and :math:`n` refers to the number of Virial parameters
:math:`a_{i}` and :math:`b_{j}`, respectively. In contrast with the isotherm models discussed earlier, Virial equation expresses the pressure
as a function of adsorption loading :math:`q^{*}` and temperature :math:`T`. Moreover, the same set of Virial parameters :math:`a_{i}` and :math:`b_{j}` are simultaneously
used to describe the isotherm data at different temperatures. The isosteric heat of adsorption :math:`Q_{st}(q^{*})` is a function of adsorption loading given as,

.. math::
  Q_{st}(q^{*}) = -R\sum_{i=0}^{m-1}a_{i}(q)^{i}

The isosteric heat of adsorption for the infinite dilution :math:`Q_{st}^{0}` (i.e. for very low adsorption loading) is expressed as,

.. math::
  Q_{st}^{0} = -R\ a_{0}

where :math:`a_{0}` is the first Virial coefficient. Note the relationship between :math:`Q_{st}` and :math:`\Delta H_{ads}` is,

.. math::
  Q_{st} = -\Delta H_{ads}

Hence, unlike Clausius-Clapeyron approach, the :math:`\Delta H_{ads}` can be directly obtained after fitting the Virial equation to
the isotherm data.

Virial Equation Fitting
~~~~~~~~~~~~~~~~~~~~~~~

In the case of Virial equation, **HeatFit** uses non-linear regression by minimizing the :math:`(SSE)` function:

.. math::
  SSE = \min_{a_{k}, b_{k}} \sum_{i}^N \left(\ln (P_{i, exp}) - f(q_{i, exp}^{*}, T_{i}; \{a_{k}\}_{0}^{m-1} \{b_{k}\}_{0}^{n-1})\right)^{2}

where :math:`N` is the total number of data points, :math:`P_{i,exp}`\, :math:`q_{i,exp}^{*}`, and :math:`T_{i}` are the experimental
pressure, experimental adsorption loading, and temperature values for the given data point :math:`i`, respectively. :math:`f` is the Virial equation; :math:`a_{k}`, :math:`b_{k}`
refers to the :math:`k^{th}` Virial coefficient, :math:`m` and :math:`n` are the total number of Virial coefficient :math:`a` and :math:`b`, respectively.

This form of :math:`SSE` is used because the Virial equation expresses the natural logarithm of adsorption pressure as a function of adsorption uptake and temperature.

**HeatFit** also reports the :math:`(RMSE)` and :math:`r^{2}` values as a measure of the goodness of fit.

.. math::
  RMSE = \sqrt{\frac{SSE}{N-M}},\\
  r^{2} = 1 - \frac{SSE}{\sum (P_{i,exp} - \overline{P_{i,exp}})}

where :math:`M` is the combined total number of Virial coefficients :math:`a` and :math:`b` and
:math:`\overline{P_{i,exp}}` is the mean value of experimental pressure.

.. note::
  Fitting multi-temperature isotherm data using Virial equation can be challenging as the same set 
  of Virial parameters :math:`a_{i}` and :math:`b_{j}` are simultaneously used to describe the isotherm data. Hence, it is important to use as few parameters as possible.
  The GUI features of **HeatFit** enable users to quickly try different numbers and combinations of Virial parameters. We advise the user to start with few parameters,
  gradually increasing the number of parameters while also keeping the check on standard errors and RMSE values, untill no significant
  improvement in RMSE is observed and the standard errors are also reasonable. Please note that parameter values with high standard error inidicate overfitting and
  are unreliable.  


.. rubric:: Footnotes
.. [#Ga] Dr. Ga Paper.