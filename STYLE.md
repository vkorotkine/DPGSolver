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
(`*_T.h` and `*_T.c`) wrapped within `def_*' and `undef_*` header files for the variable of the
specified type. For example, to define the 'd'ouble matrix functions one has in the source file:
```c
#include "def_templates_type_d.h" // type: 'd'ouble
#include "matrix_T.c"
#include "undef_templates_type.h"
```
and in the header file:
```c
#include "def_templates_type_d.h"
#include "matrix_T.h"
#include "undef_templates_type.h"
```
respectively. __Note that def/undef headers specific to the file should be included within the
templated file__ with def header files placed at the top of the file (before all declarations) and
undef header files placed at the end of the file (after all function definitions).

While this greatly reduces redundant code, it results in potential confusion when specialized
versions of the functions are called in the code despite their declarations not explicitly appearing
in any header file. Thus, it is recommended to search the documentation for function names *without*
the trailing specialization indicator (i.e. function_name`_*`) when the full function name is not
being found using `grep` for example.

Finally, in order to avoid having to provide templated names for static functions, it is recommended
to create a separate `.c` and `.h` file for each data type and simply include all of the header
files in the default `.c` and `.h` file whenever many static functions are present.

### Directory Structure

In the interest of separating unrelated ideas to the extent possible, the code has been partitioned
into several sub-directories (with `src` being the root). Due to a large number of cyclic
dependencies, many of the sub-directories with the `simulation` parent directory are linked into a
single library.

This *sometimes* complicates the build process when adding new files as it is possible to forget to
link all of the appropriate libraries or that it not be obvious which libraries may need to be
linked. In the event that an "Undefined symbols for architecture" error occurs during linking, it
may helpful to inspect the list of which functions are included as part of a static library. This
can be done by using the `nm` function in the terminal:
```sh
BUILD/path_to_lib$ nm libLIB_NAME.a
```
