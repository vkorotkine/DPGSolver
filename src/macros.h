/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__macros_h__INCLUDED
#define DPG__macros_h__INCLUDED
/** \file
 *  \brief Defines macros.
 *
 *  The call to abort() in EXIT_MSG allows for a stacktrace to be obtained when running the code using memcheck in
 *  valgrind.
 */

#include <stdio.h>
#include <stdlib.h>

///\{ \name Print the file name, function and line number.
#define PRINT_FILELINE   ({ printf("\n\nPrinting at: FILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); })
///\}

///\{ \name Exit from the code.
#define EXIT_MSG         ({ PRINT_FILELINE; fflush(stdout); abort(); })
#define EXIT_UNSUPPORTED ({printf("\n\nError: Unsupported.\n"), EXIT_MSG; })
#define EXIT_ADD_SUPPORT ({printf("\n\nError: Add support.\n"), EXIT_MSG; })
#define EXIT_DEPRECATED  ({printf("\n\nError: Deprecated.\n"), EXIT_MSG; })
#define EXIT_ERROR(...)  ({printf("\n\nError: "); printf(__VA_ARGS__); EXIT_MSG; })
///\}

///\{ \name Mark unused variables.
#define UNUSED(x)       (void)(x)
#define MAYBE_UNUSED(x) (void)(x)
///\}

///\{ \name Macro to output on successful termination.
#define OUTPUT_SUCCESS ({printf("Successful termination.\n");})
///\}

#endif // DPG__macros_h__INCLUDED
