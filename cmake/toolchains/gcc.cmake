# Toolchain file for the GCC (default)
#
# Set `CMAKE_TOOLCHAIN_FILE` while running `cmake` to use this file.


# Compilers
set(CMAKE_C_COMPILER "gcc")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED on)


# Compile flags
set(CMAKE_C_FLAGS_DEBUG          "-fPIC -g -Wall -Wextra -Werror"     CACHE STRING "Debug flags")
set(CMAKE_C_FLAGS_RELEASE        "-fPIC -O3 -DNDEBUG"                 CACHE STRING "Release flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-fPIC -O2 -g -Wall -Wextra -Werror" CACHE STRING "Release with debug info flags")
