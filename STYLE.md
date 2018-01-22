# Coding Style Guidelines

The majority of the coding style choices were selected based on the
[Linux kernel coding style guide](https://www.kernel.org/doc/html/v4.10/process/coding-style.html).

### Additional guidelines

#### General

- Tabs are NOT expanded to spaces.

- The limit on the length of lines is 120 columns.
	- 120 seemed justifiable based on the size of current terminal screens.
	- Longer lines are sometimes (but rarely) used when declaring/passing arguments to
	  functions.

- Names of user-defined types (structs) should be capitalized to distinguish them from standard
  types.
- Use conventional 'const' notation (to avoid confusion for those used to the conventional style):
```c
const int i = 0; // ok
int const i = 0; // avoid

const int*const p = NULL; // ok
int const*const p = NULL; // avoid
```

#### Header files
- What they contain:
	- Function/struct declarations;
	- The primary documentation;
	- Only **necessary** includes for the declarations.
- Nesting of header files (including headers within headers) should be avoided where possible as
  this results in a complicated dependency structure resulting in longer compile times and
  readability difficulty in assessing function interaction. Note that dependencies are however
  properly handled by CMake for compiling (i.e. CMake *will* search for implicit includes).
- This [discussion](http://stackoverflow.com/questions/1804486/should-i-use-include-in-headers)
  motivates these recommendations.

### Templating

Several functions in the code are templated such that they may be used for multiple data types
(integer, floating-point, complex floating-point), eliminating significant code duplication. In
these cases, the various function declarations/definitions are set by including the templated files
(`*_T.h` and `*_T.c`) with specialized definition (`def_*`) and undefinition (`undef_*`) macro files
in the associated header/implementation files.

While this greatly reduces redundant code, it results in potential confusion when specialized
versions of the functions are called in the code despite their declarations not explicitly appearing
in any header file. Thus, it is recommended to search the documentation for function names *without*
the trailing specialization indicator (i.e. function_name`_*`).

### Directory Structure

In the interest of separating unrelated ideas to the extent possible, the code has been partitioned into several
sub-directories (with `src` being the root). In each of the sub-directories, a static library of related functions is
generated and only the necessary libraries are then linked to the executables during the build.

This *sometimes* complicates the build process when adding new files as it is possible to forget to link all of the
appropriate libraries or that it not be obvious which libraries may need to be linked. In the event that an "Undefined
symbols for architecture" error occurs during linking, it is sometimes helpful to inspect the list of which functions
are included as part of a static library. This can be done by using the `nm` function in the terminal:
```sh
BUILD/path_to_lib$ nm libLIB_NAME.a 
```
