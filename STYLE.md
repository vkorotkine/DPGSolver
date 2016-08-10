# Coding style guide

The majority of the coding style choices were selected based on the [Linux kernel coding style
guide](https://www.kernel.org/doc/Documentation/CodingStyle)


### Additional guidelines

Header files
- What they contain: Function prototypes and minimum required headers to needed by these function prototypes.
- To ensure that the required headers have been included, the header file associated with the source file must be the
first one listed.
- Any other headers on which the source is dependent are to be included in the source file.
  - Special note: If a user-defined header is already included in the source file as part of the first header, it is
  recommended to provide a redundant include statement in the source file as well. This is done to simplify the
  dependency generation. (ToBeModified)
- This rule is broken only for main.h, which is never included in any other files.
- This [discussion](http://stackoverflow.com/questions/1804486/should-i-use-include-in-headers) motivates these
recommendations.
- To find struct header (S_*.h) dependencies hidden through the inclusion of other header files, use the following
example search command in the src directory:
  - Check if 'S_ELEMENT' is found while 'S_ELEMENT.h' is not found:
  - $ find . -name "*.c" -exec grep -lR 'S_ELEMENT' {} \; | xargs grep -L 'S_ELEMENT.h'


### Major exceptions

Tabs are 4 characters
- The only case where I have found tabs of 8 characters to be more clear than those of 4 characters is in the case of
conditional arguments which span multiple lines. Otherwise, it was felt that the longer tabs unnecessarily ate too much
into the line.

The limit on the length of lines is 120 columns
- 120 seemed justifiable based on the size of current terminal screens.
- Longer lines are sometimes used when declaring/passing arguments to functions.


Multiple statements are placed on the same line only in the following situations:
- In cases where multiple conditional/loop expressions follow each other and using a single line improves clarity;
  - While this is subjective, it seems to rarely to be the case.
- An error message is printed followed by a call to exit(1);
- Multiple free() statements are called sequentially;
- Strings are sequentially converted to long using strtol and the string pointer is moved to the end of the input;
- A linked list is searched for a specific condition;
