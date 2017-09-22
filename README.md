# Discontinuous Petrov Galerkin Solver

## Contributors

Philip Zwanenburg, philip.zwanenburg@mail.mcgill.ca

Siva Nadarajah, siva.nadarajah@mcgill.ca

Cem Gormezano

Manmeet Bhabra

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



# License

The code is licensed under the [GNU GPLv3](LICENSE.md).
