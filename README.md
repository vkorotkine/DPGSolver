# Discontinuous Petrov Galerkin Solver

### Code Description
- Uses only free to use/open source libraries/supporting programs.
- Methods:
	- Discontinuous Petrov Galerkin (DPG) (TO BE DONE);
	- Discontinuous Galerkin (DG).
- Supported elements: TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

### Supported PDEs
| PDE           | STATUS               | BRANCH (if not master) |
|---------------|----------------------| ---------------------- |
| Poisson       | DONE                 |                        |
| Advection     | DONE (NEEDS CLEANUP) | lsfem                  |
| Euler         | DONE                 |                        |
| Navier-Stokes | DONE                 |                        |


### Specific Test Cases
| Test case        |      |
|------------------|------|
| PeriodicVortex   | DONE |
| SupersonicVortex | DONE |
| InviscidChannel  | DONE |
| Taylor-Couette   | DONE |

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
| Functionality  |            |
|----------------|------------|
| MPI            | TO BE DONE |
| h/p Adaptation | DONE       |

#### Supported Numerical Fluxes
| Flux           |      |
|----------------|----- |
| Lax-Friedrichs | DONE |
| Roe-Pike       | DONE |
| Bassi-Rebay 2  | DONE |
| Compact DG 2   | DONE |

#### Supported Boundary Conditions
| Boundary Condition      |      |
|-------------------------|------|
| Dirichlet               | DONE |
| Neumann                 | DONE |
| Riemann                 | DONE |
| SlipWall                | DONE |
| BackPressure            | DONE |
| Supersonic (In/Out)flow | DONE |
| Total Temp/Pressure     | DONE |
| No Slip Dirichlet       | DONE |
| No Slip Adiabatic       | DONE |


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
