# Discontinuous Petrov Galerkin Solver

<table>
	<tr>
		<th>Branch</th>
		<th>Version</th>
		<th>Documentation</th>
	</tr>
	<tr>
		<th>
			<a href="https://github.com/PhilipZwanenburg/DPGSolver/tree/master">
				master
			</a>
		</th>
		<th>
			<a href="https://badge.fury.io/">
				<img src="https://badge.fury.io/gh/PhilipZwanenburg%2FDPGSolver.svg"
				     title="version">
			</a>
		</th>
		<th>
			<a href="https://codedocs.xyz/PhilipZwanenburg/DPGSolver/">
				<img src="https://codedocs.xyz/PhilipZwanenburg/DPGSolver.svg"
				     title="documentation">
			</a>
		</th>
	</tr>
</table>


## Code Description
- Methods:
	- Discontinuous Galerkin (DG);
	- Discontinuous Petrov Galerkin (DPG).
- Supported Partial Differential Equations: Advection, Diffusion, Euler, Navier-Stokes.
- Supported elements: LINEs, TRIs, QUADs, TETs, HEXs, WEDGEs, PYRs.
- Supported refinements: isotropic h (size) or p (order).

Note: An older version of the code can be found in the [deprecated directory](deprecated) which may
still have features not yet re-implemented in the new version, notably the 3D hp-adaptation
verification. The list of features that remains to be re-implemented can be found
[here](todo_reimplementation.md).

**It is recommended** to follow the [Coding Style Guidelines](STYLE.md) when making modifications to
the code.

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

Three separate build directories are generated based on the supported problem dimensions in which
the code can be run.

### Compile using Make

Once the code has been successfully built, the available `make` targets can be seen using
```sh
BUILD$ make help
```

Of primary interest are the following:
```sh
BUILD$ make -j     // Compile the code.
BUILD$ make meshes // Generate the meshes.
BUILD$ make -j doc // Generate the Doxygen documentation.
```

The html documentation can be accessed by pointing a browser at `BUILD/doc/html/index.html`.

### Running the Code

Executable files running various configurations of the code are placed in `BUILD/bin`, and should be executed following
the example in `$BUILD/script_files/quick.sh`. To run an executable with valgrind's memory leak detector enabled:
```sh
BUILD/script_files$ ./memcheck.sh
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

## Contributors

- Manmeet Bhabra
- Cem Gormezano
- Siva Nadarajah, siva.nadarajah (at) mcgill.ca
- Philip Zwanenburg, philip.zwanenburg (at) mail.mcgill.ca

If you would like to make your own contributions to the project, the best place to start is the 
[getting started page](https://codedocs.xyz/PhilipZwanenburg/DPGSolver/md_doc_GETTING_STARTED.html).

# License

The code is licensed under the [GNU GPLv3](LICENSE.md) due to the dependence on the GNU GSL. If this
were to be removed, the code could be relicensed under the
[MIT license](https://opensource.org/licenses/MIT).
