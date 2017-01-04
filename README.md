# Discontinuous Petrov Galerkin Solver

### Code Description
- Uses only open source libraries/supporting programs.
- Methods: Discontinuous Petrov Galerkin (DPG) with the option for standard Discontinuous Galerkin (DG).
- Supported elements: TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

### Supported PDEs
- Poisson       : DONE
- Euler         : DONE
- Navier-Stokes : TO BE DONE

See the CODE STATUS section below for details regarding current functionality.

### Installation / Set up
Follow the [installation instructions](INSTALL.md) for the set up of required libraries/programs (Currently running on
either OSX or Ubuntu). Required:
- PETSc
- MPI (MPICH or Open MPI)
- Intel (M)ath (K)ernel (L)ibrary
- ParMETIS
- Gmsh
- Paraview

Follow the instructions in [SETUP](SETUP.md) regarding additional requirements for running the code:
- Mesh generation

Please follow the [style guidelines](STYLE.md) when making additions to the code.


## Code Status
| Functionality  |            |
|----------------|------------|
| MPI            | TO BE DONE |
| h/p Adaptation | DONE       |

### Preprocessing : ACTIVE (AS NEEDED)
| Function          |                    |
|-------------------|--------------------|
| set up parameters | DONE               |
| set up mesh       | DONE               |
| set up operators  | ACTIVE (AS NEEDED) |
| set up structures | DONE               |
| set up geometry   | DONE               |

### Solving : ACTIVE
#### Initialization
| Test case        |            |
|------------------|------------|
| dSphericalBump   | TO BE DONE |
| GaussianBump     | TO BE DONE |
| PeriodicVortex   | DONE       |
| PolynomialBump   | TO BE DONE |
| SupersonicVortex | DONE       |

#### Explicit
| Function        |             |                         |
|-----------------|-------------|-------------------------|
| solver explicit |             | DONE FOR THE TIME BEING |
| volume info     |             | DONE FOR THE TIME BEING |
|                 | Weak Form   | DONE                    |
|                 | Strong Form | UNDER CONSIDERATION     |
|                 | Vectorized  | DONE                    |
| face info       |             | DONE FOR THE TIME BEING |
|                 | Weak Form   | DONE                    |
|                 | Strong Form | UNDER CONSIDERATION     |
|                 | Vectorized  | UNDER CONSIDERATION     |
| finalize        |             | DONE                    |

#### Implicit
| Function        |             |                         |
|-----------------|-------------|-------------------------|
| solver Poisson  |             | DONE                    |
| solver implicit |             | DONE FOR THE TIME BEING |
| volume info     |             | DONE FOR THE TIME BEING |
|                 | Weak Form   | DONE                    |
|                 | Strong Form | UNDER CONSIDERATION     |
|                 | Vectorized  | UNDER CONSIDERATION     |
| face info       |             | DONE FOR THE TIME BEING |
|                 | Weak Form   | DONE                    |
|                 | Strong Form | UNDER CONSIDERATION     |
|                 | Vectorized  | UNDER CONSIDERATION     |
| finalize        |             | DONE                    |

#### Fluxes
|           | Function         |                         |
|-----------|------------------|-------------------------|
| Standard  |                  | DONE FOR THE TIME BEING |
|           | inviscid         | DONE                    |
|           | viscous          | TO BE DONE              |
| Numerical |                  | DONE FOR THE TIME BEING |
|           | internal penalty | DONE                    |
|           | Bassi-Rebay 2    | DONE                    |
|           | compact DG 2     | DONE                    |
|           | Lax-Friedrichs   | DONE                    |
|           | Roe-Pike         | DONE                    |

#### Boundary
| Function  |      |
|-----------|------|
| Dirichlet | DONE |
| Neumann   | DONE |
| Riemann   | DONE |
| SlipWall  | DONE |

#### Jacobians
| Function            |      |
|---------------------|------|
| flux inviscid       | DONE |
| flux Lax-Friedrichs | DONE |
| flux Roe-Pike       | DONE |
| boundary Riemann    | DONE |
| boundary SlipWall   | DONE |


### Postprocessing : ACTIVE (As cases are added)


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
