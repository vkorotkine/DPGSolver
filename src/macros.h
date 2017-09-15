// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__macros_h__INCLUDED
#define DPG__macros_h__INCLUDED
/** \file
 *	\brief Defines macros.
 *
 *	The call to abort() in EXIT_MSG allows for a stacktrace to be obtained when running the code using memcheck in
 *	valgrind.
 */

#include <stdio.h>
#include <stdlib.h>

///\{ \name Print the file name, function and line number.
#define PRINT_FILELINE   ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); })
///\}

///\{ \name Exit from the code.
#define EXIT_MSG         ({ PRINT_FILELINE; abort(); })
#define EXIT_UNSUPPORTED ({printf("\n\nError: Unsupported.\n"), EXIT_MSG; })
#define EXIT_ADD_SUPPORT ({printf("\n\nError: Add support.\n"), EXIT_MSG; })
#define EXIT_ERROR(...)  ({printf("\n\nError: "); printf(__VA_ARGS__); EXIT_MSG; })
///\}

///\{ \name Mark unused variables.
#define UNUSED(x) (void)(x)
///\}

#endif // DPG__macros_h__INCLUDED
