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
