# Discontinuous Petrov Galerkin Solver

### Code Description
- Uses only free to use/open source libraries/supporting programs.
- Methods:
	- Discontinuous Galerkin (DG).
	- Discontinuous Petrov Galerkin (DPG) (TO BE DONE);
- Supported elements: TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

### Supported PDEs
| PDE           | STATUS               | COMMENTS |
|---------------|----------------------|----------|
| Advection     | DONE                 |Additional functionality in lsfem branch|
| Poisson       | DONE                 ||
| Euler         | DONE                 ||
| Navier-Stokes | DONE                 ||


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

See the CODE STATUS section below for details regarding current functionality.

### Installation / Set up
Follow the [installation instructions](INSTALL.md) for the set up of required libraries/programs. Required:
- MPI (MPICH or Open MPI)
- Intel (M)ath (K)ernel (L)ibrary
- PETSc
- ParMETIS
- Gmsh
- Paraview
- Python3

Follow the instructions in [SETUP](SETUP.md) regarding additional requirements for running the code:
- Mesh generation

Please follow the [style guidelines](STYLE.md) when making additions to the code.


## Code Status
| Functionality  | Status     |
|----------------|------------|
| MPI            | TO BE DONE |
| h/p Adaptation | DONE       |

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


### License
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
