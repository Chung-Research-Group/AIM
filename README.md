<p align="center">
  <img src="src/AIM_logo.png" alt="Logo" width="300"/>
</p>

AIM
========

AIM (Adsorption Integrated Modules) is a collection of MATLAB based Graphical User Interface (GUI) tools for adsorption isotherm based fixed bed modelling. AIM modules provide an integrated workflow for adsorpion isoterm fitting, isosteric heat of adsorption estimation, mixture isotherm prediction, and multicomponent adsorption breakthrough simulation. AIM has been developed at [Chung Group](https://sites.google.com/view/mtap-lab), Pusan National University.

Primary features of AIM include
- Integrated workflow within modules
- GUI environnment to facilitate data input/output, model selection, and simulation type selection
- Dynamic plots for pressure, loading, temperature, and composition
- Data and graphics export in different formats

**Developed by:** [Muhammad Hassan](https://github.com/hassan-azizi)

Installation
========
- **Windows users**:
   <br>Please use the 'AIM_Installer_v_1.0.1.exe' file located in the [bin](./bin/) directory.

- **Linux and macOS users**:
   <br>Users of Linux and macOS need to compile the app for redistirbution.
   <br>Please follow the instruction given in the [README](./build/README.md) file located in the build directory for compilation.

AIM Modules
======
- **IsoFit**: Single temperature isotherm fitting using various isotherm models
   * Single and dual-site Langmuir
   * Single and dual-site Langmuir-Freundlich
   * Quadratic
   * Temkin's approximation
   * BET
   * Sips
   * Toth
   * Structural Transition Adsorption
   * Dubinin-Astakhov
   * Klotz
   * Do-Do

- **HeatFit**: Multiple temperature isotherm fitting and isosteric heat of adsorption prediction using
   * Clausius-Clapeyron
   * Virial equation

- **MixPred**: Mixture adsorption isotherm using
   * Extended Dual-site Langmuir model (EDSL)
   * Ideal Adsorption Solution Theory (IAST)

- **BreakLab**: Rigorous multicomponent breakthrough simulation
   * Supports up to 5 component systems
   * Non-isothermal/non-isobaric breaktrhough simulation
   * Axial dispersion
   * Linear Driving Force (LDF)
   * Ergun equation for pressure drop

Cite Us
============
If you use AIM software for your scientific publications, please cite:<br>
**"AIM: A User-friendly GUI Workflow program for Isotherm Fitting, Mixture Prediction, Isosteric Heat of Adsorption Estimation, and Breakthrough Simulation"**<br>
Muhammad Hassan, Sunghyun Yoon, Yu Chen, Pilseok Kim, Hongryeol Yun, Youn-Sang Bae, Chung-Yul Yoo, Dong-Yeun Koh, Chang-Seop Hong, Ki-Bong Lee, Yongchul G. Chung<br>
Journal: ####<br>
URL: ####

# License
This project is licensed under the GNU General Public License v2 (GPLv2).  
See the [LICENSE](./LICENSE) file for details.

# Maintainers
* Muhammad Hassan (hassanaz.14@outlook.com)
