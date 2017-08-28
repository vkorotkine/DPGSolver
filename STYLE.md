# Coding Style Guidelines

The majority of the coding style choices were selected based on the [Linux kernel coding style guide](https://www.kernel.org/doc/html/v4.10/process/coding-style.html).

### Additional guidelines

#### General

- Tabs are NOT expanded to spaces.

- The limit on the length of lines is 120 columns.
	- 120 seemed justifiable based on the size of current terminal screens.
	- Longer lines are sometimes (but rarely) used when declaring/passing arguments to functions.

##### May not currently be the case ... Correct where occurences are found
- Names of user-defined types (structs) should be capitalized to distinguish them from standard types.
- Use conventional 'const' notation (to avoid confusion for those used to the conventional style):
```c
const int i = 0; // ok
int const i = 0; // avoid

const int*const p = NULL; // ok
int const*const p = NULL; // avoid
```



#### Header files
- What they contain: Function/struct declarations and constant definition headers needed by these declarations.
- Nesting of header files should be avoided where possible (including headers within headers) as this results in
  complicated dependency treatment which is difficult to track correctly using the current dependency generation
  mechanism. If a header is needed in a file (.c) and is included indirectly through a header, it will **not** be
  present in the dependency file (.d).
	- I believe that the nested dependencies can be handled by CMake which should be used in future, although I have not
	  tested this extensively.
	- This results in:
		- repeated declaration of structs, but limits hidden dependency propagation.
	- The *only* exceptions to this rule are:
		- the inclusion of headers defining constants (as they cannot be declared extern). In this case, any constant
		  definition headers should be redundantly reincluded in the associated '.c' file immediately following the
		  primary header inclusion.
		- the inclusion of header files defining structs which are members of structs declared in the header. **Ensure
		  that the nested dependencies are redundantly included in all header (.h) and source (.c) files!**.
- This [discussion](http://stackoverflow.com/questions/1804486/should-i-use-include-in-headers) motivates these
recommendations.
