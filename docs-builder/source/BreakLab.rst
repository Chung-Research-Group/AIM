.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BreakLab
===============================

**BreakLab** is the module for isothermal/non-isothermal multicomponent fixed bed breakthrough simulation for up to 5 components.
**BreakLab** GUI  is shown below:  
 
.. image:: images/BreakLab.png
   :width: 1000
   :alt: BreakLab_logo
   :align: center
   
-------

The isotherm fitting results from **IsoFit** and **HeatFit** for different components can be directly loaded into **BreakLab** module.

Watch how to use **BreakLab** :ref:`here <BreakLab-label>`

The mathematical model implemented in BreakLab and the required properties and parameters for breakthrough simulation are summarised below.

Mathematical Model
---------------------------------------------

The mathematical model of BreakLab is based on the following assumptions:

* The gas flow is axially dispersed and charecterized by an axial dispersion coefficient.
* The gas phase is ideal.
* Pressure drop in the column is given by Ergun's equation.
* The steady state momentum balance is applicable.
* The adsorbent bed is uniform; the bed density, void fraction, and particle size are constant throughout the bed length.
* Thermal equilibrium exists between gas and solid phase.
* The heat transfer coefficient governing the heat transfer between column wall and bed is constant.
* The wall temperature remains constant.
* The mass transfer resistance between solid and gas phase is governed by Linear Driving Force (LDF) model.
* The concentration, pressure, and temperature gradients in the axial directions are negligible.

The material, energy and momentum balances are developed based on conservation of mass, energy, and momentum in the fixed bed, respectively.
The resulting balance equations are partial differential equations (PDEs) and consist of mole conservation of component :math:`i`, overall mass balance, conservation of total momentum and the total energy
of the system. The balance equations are summarised in table below. Please check the journal article for a more detailed discussion on deriving these 
forms of balance equations.

.. list-table:: BreakLab Balance Equations
   :header-rows: 1
   :widths: auto

   * - 
     - Expression
   * - Component mole balance
     - :math:`\varepsilon_{t}\left(\frac{\partial y_{i}}{\partial t} + \frac{y_{i}}{P}\frac{\partial P}{\partial t} - \frac{y_{i}}{T}\frac{\partial T}{\partial t} \right) = D_{ax} \frac{\varepsilon_{b} T}{P}\frac{\partial}{\partial z} \left(\frac{P}{T} \frac{\partial y_{i}}{\partial z}\right) - \frac{\varepsilon_{b} T}{P}\frac{\partial}{\partial z} \left(\frac{y_{i}vP}{T}\right) - \frac{\rho_{b,ads}RT}{P} \frac{\partial q_{i}}{\partial t}`
   * - Total mole balance
     - :math:`\frac{\partial P}{\partial t} = \frac{P}{T}\frac{\partial T}{\partial t} - \frac{\varepsilon_{b} T}{\varepsilon_{t}}\frac{\partial}{\partial z} \left(\frac{vP}{T}\right) - \frac{\rho_{b,ads}RT}{\varepsilon_{t}} \sum_{i\in I} \frac{\partial q_{i}}{\partial t}`
   * - Energy Balance
     - :math:`\left(\rho_{b,ads}C_{p,ads} + \rho_{b,ads}C_{p,a}\sum_{i \in I} q_{i}\right)\frac{\partial T}{\partial t} = K_{z} \frac{\partial^{2}T}{\partial z^2} - \frac{C_{p,gas}\varepsilon_{b}}{R}\frac{\partial}{\partial z}(vP) + \rho_{b,ads} \left( \sum_{i \in I} \left(-\Delta H_{ads,i}\frac{\partial q_{i}}{\partial t}\right) - C_{p, a}T \sum_{i \in I}\frac{\partial q_{i}}{\partial t}\right) - \frac{2 h_{in}}{r_{in}} (T - T_{wall}) - \frac{C_{p,gas}\varepsilon_{t}}{R}\frac{\partial P}{\partial t}`
   * - Linear Driving Force
     - :math:`\frac{\partial q}{\partial t} = k_{i}(q_{mix, i}^{*} - q_{i})`
   * - Ergun Equations
     - :math:`- \frac{\partial P}{\partial z} = \left(\frac{150 \mu}{4r_{p}^{2}}\right) \left(\frac{1-\varepsilon_{b}}{\varepsilon_{b}}\right)^2v + \left(\frac{1.75 \rho_{gas}}{2r_{p}} \left(\frac{1-\varepsilon_{b}}{\varepsilon_{b}}\right) \right)v^{2}` 


Total porosity :math:`\varepsilon_{t}` is calculated using the bulk and paricle porosities, :math:`\varepsilon_{b}`, :math:`\varepsilon_{p}`,

.. math::
  \varepsilon_{t} =  \varepsilon_{b} + (1-\varepsilon_{b})\varepsilon_{p}

Axial dispersion coefficient :math:`D_{ax}` is calculated using,

.. math::
  D_{ax} = 0.7D_{m} + v_{0}r_{p}

Boundary Conditions
~~~~~~~~~~~~~~~~~~~~

BreakLab uses the Danckwerts boundary conditions (BCs) for dispersed plug flow system. The BCs are summarised in table below.

.. list-table:: BreakLab Boundary Conditions
   :header-rows: 1
   :widths: auto

   * - 
     - Inlet (:math:`z = 0`)
     - Outlet (:math:`z = L`)
   * - Mole fraction (:math:`y_{i}`)
     - :math:`D_{ax} \frac{\partial y_{i}}{\partial z}|_{z=0} = -v|_{z=0}(y_{0, i} - y_{i}|_{z=0})`
     - :math:`\frac{\partial y_{i}}{\partial z}|_{z=L} = 0`
   * - Temperature (:math:`T`)
     - :math:`K_{z} \frac{\partial T}{\partial z}|_{z=0} = -\varepsilon_{b}\rho_{gas}C_{p, gas}v|_{z=0}(T_{0} - T|_{z=0})`
     - :math:`\frac{\partial T}{\partial z}|_{z=L} = 0`
   * - Velocity (:math:`v`)
     - :math:`v|_{z=0} = v_{0}/\varepsilon_{b}`
     - None
   * - Pressure (:math:`P`)
     - :math:`P|_{z=0} = f^*(v|_{z=0})`
     - :math:`P|_{z=L} = P_{0}`

\* here :math:`f` represents the Ergun equation. Only the velocity is specified at inlet and pressure is back calculated using Ergun equation.

Together the balance equations and BCs consitute the mathematical model of BreakLab.

Numerical Simulation
---------------------------------------------

The model equations are nondimensionalized to resolve the steep gradients in the adsorption processes and faster convergence. The scaling variables
used for nondimensionalization are as follows:

.. math::
  \overline{P} = \frac{P}{P_{0}},\quad \overline{T} =\frac{T}{T_{0}},\quad \overline{v} = \frac{v}{v_{0}}, \quad \overline{x_{i}} = \frac{q_{i}}{q_{i0}},
  \quad {\xi} = \frac{z}{L},\quad \tau = \frac{tv_{0}}{L}, \quad \overline{T_{wall}} = \frac{T_{wall}}{T_{0}}

The dimensionless equations are derived using the scaled variables above. BreakLab employs finite volume method (FVM) along with weighted essentially non-oscillatory (WENO) scheme
for spatial discretization of the above PDEs. Spatial discretization yields ordinary differential equations (ODEs). BreakLab uses MATLAB's built in stiff ODE solver *ode15s* for solving the systems of ODEs.

Breakthrough Simulation Input
-----------------------------

The table below summarises all the input properties and parameters required for running a breakthrough simulation in BreakLab.

.. list-table:: BreakLab Inputs
   :header-rows: 1
   :widths: auto

   * - Parameters
     - Description
     - Units
     - Remarks
   * - :math:`T_{0}`
     - Feed gas temperature
     - :math:`^\circ\mathrm{C}`
     - :math:`-`
   * - :math:`v_{0}`
     - Superficial feed gas velocity
     - :math:`m/s\\\mathrm{OR}\\\mathrm{std. }\quad cm^3/min`
     - Standard state refers to :math:`0\,^\circ\mathrm{C}` and :math:`101325\,Pa`
   * - :math:`y_{0, i}`
     - Feed gas mole fraction of component :math:`i`
     - :math:`-`
     - :math:`-`    
   * - :math:`MW_{i}`
     - Molecular weight of component :math:`i`
     - :math:`kg\,mol^{-1}`
     - :math:`-`
   * - :math:`D_{m}`
     - Molecular diffusivity of gas
     - :math:`m^{2}/s`
     - :math:`-`
   * - :math:`K_{z}`
     - Thermal conductivity of gas
     - :math:`W/m/^\circ\mathrm{C}`
     - :math:`-`
   * - :math:`C_{p, gas}`
     - Specific heat capacity of gas
     - :math:`J/mol/^\circ\mathrm{C}`
     - :math:`-`
   * - :math:`\mu`
     - Viscosity of gas
     - :math:`Pa\,s`
     - :math:`-`
   
   * - :math:`P_{0}`
     - Absolute pressure of adsorber column
     - :math:`kPa`
     - This pressure refers to the pressure at the outlet of column.
   * - :math:`L`
     - Length of adsorber column
     - :math:`m`
     - :math:`-`
   * - :math:`D_{col}`
     - Inner diameter adsorber column
     - :math:`m`
     - :math:`-`
   * - :math:`h_{in}`
     - Inside heat transfer coefficient
     - :math:`W/m^{2}/^\circ\mathrm{C}`
     - For running adiabatic breakthrough simulations specify this equal to :math:`0.0`.
   * - :math:`T_{wall}`
     - Temperature of adsorber wall
     - :math:`^\circ\mathrm{C}`
     - :math:`-`

   * - :math:`\rho_{b,ads}`
     - Bulk density of adsobent
     - :math:`kg/m^3`
     - :math:`-`
   * - :math:`D_{p}`
     - Diameter of adsorbent particles
     - :math:`m`
     - :math:`-`
   * - :math:`\varepsilon_{b}`
     - Bulk porosity of adsorbent bed
     - :math:`-`
     - :math:`-`
   * - :math:`\varepsilon_{p}`
     - Particle porosity of adsorbent particles
     - :math:`-`
     - In case of not knowing the particle porosity specify this as :math:`0.0`.
   * - :math:`C_{p,ads}`
     - Specific heat capacity of adsorbent
     - :math:`J/kg/^\circ\mathrm{C}`
     - :math:`-`
   
   * - Isotherm parameters
     - Isotherm parameters for the chosen isotherm model
     - :math:`-`
     - Need to provide this for each component except for carrier gas. \
       
       The number and types of isotherm model parameters depends on the 
       
       chosen isotherm model. The isotherm model and parameters for the \
       
       given component obtained from IsoFit or HeatFit can be imported as well. 
   * - :math:`\Delta H_{ads, i}`
     - Heat of adsorption of component :math:`i`
     - :math:`kJ/mol`
     - :math:`-`
   * - :math:`T_{ref}`
     - Reference temperature used for isotherm fitting
     - :math:`K`
     - :math:`-`
   * - :math:`k_{i}`
     - Mass transfer coefficient of component :math:`i`
     - :math:`1/s`
     - :math:`-`
   
            
Nomenclature
-------------

| :math:`P =` Absolute pressure of the gas :math:`(Pa)`
| :math:`\overline{P} =` Dimensionless pressure of the gas :math:`(-)`
| :math:`T =` Absolute Temperature of the gas :math:`(K)`
| :math:`\overline{T} =` Dimensionless temperature of the gas :math:`(-)`
| :math:`T_{wall} =` Wall temperature :math:`(K)`
| :math:`\overline{T_{wall}} =` Dimensionless wall temperature :math:`(-)`
| :math:`y_{i} =` Gas phase mole fraction of component :math:`i` :math:`(Pa)`
| :math:`q_{i} =` Adsorption loading of component :math:`i` :math:`(mol/kg)`
| :math:`\overline{x_{i}} =` Dimensionless adsorption loading of component :math:`i` :math:`(-)`
| :math:`q_{i}^{*} =` Equilibrium adsorption loading of component :math:`i`  :math:`(mol/kg)`
| :math:`L =` Adsorber column length :math:`(m)`
| :math:`\xi =` Dimensionless adsorber column length :math:`(-)`
| :math:`\overline{t} =` Time :math:`(s)`
| :math:`\tau =` Dimensionless time :math:`(-)`
| :math:`\varepsilon_{t} =` Total porosity of adsorbent particle :math:`(-)`
| :math:`\varepsilon_{b} =` Bulk porosity of adsorbent particle :math:`(-)`
| :math:`\varepsilon_{p} =` Particle porosity of adsorbent particle :math:`(-)`
| :math:`v =` Interstitial velocity :math:`(m/s)`
| :math:`v_{0} =` Superficial velocity of feed :math:`(m/s)`
| :math:`P_{0} =` Adsorption pressure :math:`(Pa)`
| :math:`T_{0} =` Temperature of feed :math:`(K)`
| :math:`y_{0, i} =` Mole fraction of component :math:`i` in feed :math:`(-)`
| :math:`\rho_{b, ads} =` Bulk density of adsorbent bed :math:`(kg/m^{3})`
| :math:`\rho_{gas} =` Density of gas phase :math:`(kg/m^{3})`
| :math:`D_{ax} =` Axial dispersion coefficient :math:`(m^{2}/s)`
| :math:`D_{m} =` Molecular diffusivity of gas phase :math:`(m^{2}/s)`
| :math:`r_{p} =` Radius of adsorbent particles :math:`(m)`
| :math:`r_{in} =` Inside radius of adsorber column :math:`(m)`
| :math:`\mu =` Viscosity of gas phase :math:`(Pa . s)`
| :math:`K_{z} =` Thermal conductivity :math:`(W/m/K)`
| :math:`C_{p,ads} =` Heat capacity of adsorbent particles :math:`(J/kg/K)`
| :math:`C_{p,a} =` Heat capacity of adsorbed phase :math:`(J/mol/K)`
| :math:`C_{p,gas} =` Heat capacity of gas phase :math:`(J/mol/K)`
| :math:`\Delta H_{ads, i} =` Heat of adsorption of component :math:`i` :math:`(J/mol)`
| :math:`h_{in} =` Inside heat transfer coefficient :math:`(W/m^{2}/K)`
| :math:`k_{i} =` Mass transfer coefficient of component :math:`i`  :math:`(1/s)`