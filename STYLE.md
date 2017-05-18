# Coding style guide

The majority of the coding style choices were selected based on the [Linux kernel coding style
guide](https://www.kernel.org/doc/html/v4.10/process/coding-style.html)

### Additional guidelines

Header files
- What they contain: Function prototypes and minimum required headers needed by these function prototypes.
- To ensure that the required headers have been included, the header file associated with the source file must be the
first one listed.
- Any other headers on which the source is dependent are to be included in the source file.
  - Special note: If a user-defined header is already included in the source file as part of the first header, it is
  recommended to provide a redundant include statement in the source file as well. This is required for the dependency
  generation.
- This [discussion](http://stackoverflow.com/questions/1804486/should-i-use-include-in-headers) motivates these
recommendations.
- To find struct header (S_*.h) dependencies hidden through the inclusion of other header files, use the following
example search command in the src directory:
  - Check if 'S_ELEMENT' is found while 'S_ELEMENT.h' is not found:
  - $ find . -name "*.c" -exec grep -lR 'S_ELEMENT' {} \; | xargs grep -L 'S_ELEMENT.h'

Multiple statements placed on the same line should generally be avoided but may be u

The limit on the length of lines is 120 columns
- 120 seemed justifiable based on the size of current terminal screens.
- Longer lines are sometimes used when declaring/passing arguments to functions.
