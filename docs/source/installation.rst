Installation
============

Windows standalone application
------------------------------

The repository contains a Windows installer at
``bin/AIM_Installer_v_1.0.exe``.

1. Download the installer from the
   `bin directory <https://github.com/Chung-Research-Group/AIM/tree/main/bin>`_.
2. Run it and follow the installation prompts.
3. Install the matching MATLAB Runtime if prompted.

The compiled application does not require a MATLAB license. The supplied
installer is Windows-only.

Run from source
---------------

Source use requires MATLAB R2024a or newer.

.. code-block:: console

   git clone https://github.com/Chung-Research-Group/AIM.git
   cd AIM

Open ``AIM_Build.prj`` in MATLAB, then run ``src/Main_app.mlapp``. Keep the
project file at the repository root so that project paths resolve correctly.

Build a standalone application
------------------------------

Packaging requires MATLAB Compiler R2024a or newer. Open ``AIM_Build.prj`` and
use MATLAB's **Package** action. Standalone builds are platform-dependent and
must use a MATLAB Runtime version matching the MATLAB Compiler release and
update used for packaging.

See the repository's
`build guide <https://github.com/Chung-Research-Group/AIM/blob/main/build/README.md>`_
for deployment details.
