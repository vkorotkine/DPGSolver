# Coding style guide

The majority of the coding style choices were selected based on the [Linux kernel coding style
guide](https://www.kernel.org/doc/Documentation/CodingStyle)

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
