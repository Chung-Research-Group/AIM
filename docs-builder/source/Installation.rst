.. AIM Documentation documentation master file, created by
   sphinx-quickstart on Fri May 16 14:38:34 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
===============================

Windows
-----------------

AIM GUI software can be freely downloaded from :download:`here <_static/AIM-Installer.tar.gz>`

1. Unzip the ``AIM-Installer.tar.gz`` file.
2. Use ``AIM_Installer.exe`` file for installing AIM.

AIM application requires MATLAB Runtime environment. The AIM Installer will automatically install both AIM application and MATLAB Runtime.

Users can also specify the desired installation location for both AIM and MATLAB Runtime.

Please note that ``AIM_Installer.exe`` file will work only on Windows platform.

Linux and macOS
-----------------

Users of Linux and macOS need to first compile AIM, if they want to use or re-distribute it as a standalone application.
Please follow the instructions below to compile the provided source code using the MATLAB for re-distribution. 
The successful compilation will generate AIM installer that can be distributed. Note that the MATLAB compilation is platform dependent. 
The generated executable will function properly only if the operating system of the end user is same as the operating system of the machine, the code is compiled on.

Prerequisites
~~~~~~~~~~~~~

**MATLAB Version**: Ensure you have MATLAB R2024a or a newer version installed.

**MATLAB Compiler**: This MATLAB toolbox is required to compile MATLAB code into standalone applications. Version must be R2024a or newer.

**MATLAB Runtime**: The compiled application requires the MATLAB Runtime environment to run on machines without MATLAB installed. The version of MATLAB Runtime must match the version of MATLAB Compiler used for compilation. 
While compiling the source code, one has two options:

1. Create an application installer with MATLAB Runtime.
2. Create an application installer without MATLAB Runtime. In this case, users of the application will need to download the MATLAB Runtime.

In case 1, you must download the MATLAB Runtime installer that matches both the version and update level of MATLAB used to create the installer.
Details about instaling MATLAB Runtime on different platforms can be found on `MATLAB Website <https://www.mathworks.com/help/compiler/install-the-matlab-runtime.html>`_

Compilation
~~~~~~~~~~~~~

1. **Download the AIM sourcecode**
   
   * The source code of AIM is opensource and is freely avaibale at our GitHub repository `AIM Repository <https://github.com/mtap-research/AIM>`_.
     You can clone the repository using Git:

      .. code-block:: bash

         git clone https://github.com/mtap-research/AIM.git

2. **Start MATLAB**
   
   Start the MATLAB application.

3. **Open AIM Project**

   * Navigate to the downloaded AIM repository directory within MATLAB.
   * Double click the ``AIM_Build.prj`` file to open it in MATLAB.

4. **Check Necessary Files**
   
   * The project template file ``AIM_Build.prj`` is already configured and contain all the neccessary files.
   * Do not change the location of the ``AIM_Build.prj`` file. It must be present in the AIM repository root directory.

5. **Compile the Application**
   
   * Click the "Package" button to start the compilation process.
   * Upon completion, MATLAB will generate the standalone application and a ``readme.txt`` file containing deployment prerequisites and a list of packaged files.

