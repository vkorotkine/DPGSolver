# Discontinuous Petrov Galerkin Solver

### Code Description
- Open source using only open source libraries/supporting programs.
- Methods: Discontinuous Petrov Galerkin (DPG) with the option for standard Discontinuous Galerkin (DG).
- Supported elements: LINEs, TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinement: anisotropic (LINE/QUAD/HEX)/ isotropic (TRI/TET/WEDGE/PYR) h (size), p (order).

### Supported PDEs
- Euler         (ACTIVE)
- Navier-Stokes TO BE DONE

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
- MPI solving not yet supported.

### Preprocessing  : ACTIVE
- set up parameters : DONE
- set up mesh       : DONE
- set up operators  : ACTIVE (Solver)
- set up structures : ACTIVE (Solver)
- set up geometry   : DONE

### Solving        : ACTIVE
#### Initialization
- dSphericalBump   : TO BE DONE
- GaussianBump     : TO BE DONE
- PeriodicVortex:  : DONE
- PolynomialBump   : TO BE DONE
- SupersonicVortex : TO BE DONE

#### Explicit
- solver RK : ACTIVE
- VOLUME    : ACTIVE
  - Weak Form   : DONE
  - Strong Form : TO BE DONE
  - Vectorized  : TO BE DONE
- FACET     : ACTIVE
  - Weak Form   : DONE
  - Strong Form : TO BE DONE
  - Vectorized  : TO BE DONE
- finalize  : DONE

#### Implicit
- solver implicit TO BE DONE
- VOLUME          TO BE DONE
- FACET           TO BE DONE
- finalize        TO BE DONE

#### Fluxes
- standard
  - inviscid       : DONE
  - viscous        : TO BE DONE
- numerical
  - Lax-Friedrichs : DONE
  - Roe-Pike       : DONE

#### Boundary
- riemann      : TO BE DONE
- slip wall    : TO BE DONE
- outflow mach : TO BE DONE

#### Jacobians
- flux inviscid         : TO BE DONE
- flux LF               : TO BE DONE
- flux Roe              : TO BE DONE
- boundary riemann      : TO BE DONE
- boundary slip wall    : TO BE DONE
- boundary outflow mach : TO BE DONE


### Postprocessing : TO BE DONE


### License
The MIT License (MIT)

Copyright (c) 2016 Philip Zwanenburg

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
