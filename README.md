# Discontinuous Petrov Galerkin Solver

## Code Description
- Uses only free to use/open source libraries/supporting programs.
- Methods:
	- Discontinuous Galerkin (DG);
	- Hybridized Discontinuous Galerkin (HDG) (ACTIVE);
	- Discontinuous Petrov Galerkin (DPG) (TO BE DONE);
- Supported elements: LINEs, TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

Please follow the [Coding Style Guidelines](STYLE.md) when making modifications to the code.

## Code Status

### General
| Functionality  | Status     |
|----------------|------------|
| MPI            | TO BE DONE |
| h/p Adaptation | DONE       |


### Test Cases
| PDE           | Name             | Status |
|---------------|------------------|--------|
| Advection     | Default          | DONE   |
|               | Peterson         | DONE   |
| Poisson       | Default          | DONE   |
| Euler         | PeriodicVortex   | DONE   |
|               | SupersonicVortex | DONE   |
|               | InviscidChannel  | DONE   |
| Navier-Stokes | Taylor-Couette   | DONE   |


#### Supported Numerical Fluxes
| Name           | Status |
|----------------|------- |
| Upwind         | DONE   |
| Lax-Friedrichs | DONE   |
| Roe-Pike       | DONE   |
| Bassi-Rebay 2  | DONE   |
| Compact DG 2   | DONE   |


#### Supported Boundary Conditions
| PDE           | Name                    | Status |
|---------------|-------------------------|--------|
| Advection     | Inflow                  | DONE   |
|               | Outflow                 | DONE   |
| Poisson       | Dirichlet               | DONE   |
|               | Neumann                 | DONE   |
| Euler         | Riemann                 | DONE   |
|               | SlipWall                | DONE   |
|               | BackPressure            | DONE   |
|               | Supersonic (In/Out)flow | DONE   |
|               | Total Temp/Pressure     | DONE   |
| Navier-Stokes | No Slip Dirichlet       | DONE   |
|               | No Slip Adiabatic       | DONE   |



## Installation / Set up
Follow the [Detailed Installation Instructions](INSTALL.md) for the set up of **required** external packages.

Follow the instructions in [Additional Set Up](SETUP.md) regarding additional requirements for running the code:
- Configuring the code
- Running the code
- Generating documentation
	- Documentation is very incomplete. \todo Add/update documentation where necessary.

#### External Packages

##### Required
- [MPICH](https://www.mpich.org)
- [Intel (M)ath (K)ernel (L)ibrary](https://software.intel.com/en-us/mkl)
- [PETSc](https://www.mcs.anl.gov/petsc/)
- [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
- [Gmsh](http://gmsh.info)
- [Python3](https://www.python.org/downloads/) (including [Numpy](http://www.numpy.org))

##### Visualization
- [Paraview](https://www.paraview.org)

##### Documentation
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/)
- [graphivz](http://www.graphviz.org)




## License
The MIT License (MIT)

Copyright (c) 2017 Philip Zwanenburg

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
