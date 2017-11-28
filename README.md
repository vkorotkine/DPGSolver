# Discontinuous Petrov Galerkin Solver

## Contributors

Philip Zwanenburg, philip.zwanenburg@mail.mcgill.ca

Siva Nadarajah, siva.nadarajah@mcgill.ca

Cem Gormezano

Manmeet Bhabra

## Building/Running the Code

The code has been succesfully built in the following environments:
- linux (ubuntu 16.04);
- macOS (Sierra 10.12.5).

### Build using CMake

An out-of-source build must be performed using the [sample scripts](cmake/run) by executing the
appropriate bash script. For example, to configure for the debug build for macOS:
```sh
$ ROOT/cmake/run$ sh macOS_gcc_debug.sh
```

**A customized script file may be required** if CMake is unable to locate required software which
you are sure is installed; CMake will exit with an error if required software is not found. It is
recommended to use a standard package manager for installing missing software whenever possible
(e.g. `apt` on ubuntu and `homebrew` on macOS) to limit package conflicts and place installed
software in directories in the default search path.

Three separate build directories are generated based on the supported problem dimensions in which the code can be run.

### Compile using Make

Once the code has been successfully built, the available `make` targets can be seen using
```sh
BUILD$ make help
```

Of primary interest are the following:
```sh
BUILD$ make        // Compile the code.
BUILD$ make meshes // Generate the meshes.
BUILD$ make doc    // Generate the Doxygen documentation.
```

The html documentation can be accessed by pointing a browser at `BUILD/doc/html/index.html`.

### Running the Code

Executable files running various configurations of the code are placed in `BUILD/bin`, and should be executed following
the example in `$BUILD/script_files/quick.sh`. To run an executable with valgrind's memory leak detector enabled:
```sh
BUILD/script_files$ sh memcheck.sh
```

## Testing

**Note: In order to run all of the available tests, the code must be compiled in all of the dimension-dependent build
directories.**

Testing can be performed using CMake's `ctest` functionality. After successfully compiling the project, all tests can be
run by executing:
```sh
BUILD$ ctest (which is equivalent to BUILD$ make test)
```

An alternative make target is provided to run tests with --output-on-failure:
```sh
BUILD$ make check
```

Additional useful commands are:
```sh
BUILD$ ctest -N (List the tests that would be run but not actually run them)
BUILD$ ctest -R <regex> (Run tests matching regular expression)
BUILD$ ctest -V (Enable verbose output from tests)
```

## Code Description
- Methods:
	- Discontinuous Galerkin (DG);
	- Hybridized Discontinuous Galerkin (HDG);
	- Discontinuous Petrov Galerkin (DPG).
- Supported elements: LINEs, TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

**Please** follow the [Coding Style Guidelines](STYLE.md) when making modifications to the code.

## Code Status

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

# License

The code is licensed under the [GNU GPLv3](LICENSE.md) due to the dependence on the GNU GSL. If this
were to be removed, the code could be relicensed under the
[MIT license](https://opensource.org/licenses/MIT).
