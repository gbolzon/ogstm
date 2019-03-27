## A few notes about the documentation

#### Template namelists

A few commented template namelists are provided as an example:

- [`./doc/boundaries.nml.template`](./doc/boundaries.nml.template): template namelist with the main list of boundaries;
- [`./doc/riv.nml.template`](./doc/riv.nml.template): template namelist for a `rivers` boundary;
- [`./doc/gib.nml.template`](./doc/gib.nml.template): template namelist for a `sponge` boundary with nudging;
- [`./doc/dar.nml.template`](./doc/dar.nml.template): template namelist for an `open` boundary.

#### Additional documentation

Further insight on the motivation and a detailed description of the boundary conditon modules can be found in the following documents:

- [`./doc/marco_bettiol_MHPC_2017-18_thesis.pdf`](./doc/marco_bettiol_MHPC_2017-18_thesis.pdf): main thesis work developed under the MHPC 2017-18 program, referencing to commit `e19e5d` ("complete documentation");
- [`./doc/marco_bettiol_MHPC_2017-18_thesis_pres_interface.pdf`](./doc/marco_bettiol_MHPC_2017-18_thesis_pres_interface.pdf): presentation of the new interface for the boundary conditions, with reference to commit `cee8bd42` (Bug fix).

#### `Doxygen`

Full documentation of the code can be auto-generated using `Doxygen`. A Doxygen configuration file is alredy provided ([`./src/doxygen.conf`](./src/doxygen.conf)) and two documentation formats (`latex` and `html`) are supported. Documentation source files in both formats are provided with `doxygen doxygen.conf`. Final documentaion can be then obtained calling `make` within the desider subfolder. Tests have been run with `Doxygen 1.8.14`. Documentation outputs are reported in [`./src/doxygen/`](./src/doxygen/) (directory may need to be created first).