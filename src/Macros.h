// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Macros_h__INCLUDED
#define DPG__Macros_h__INCLUDED

/*
 *	Purpose:
 *		Define macros.
 *
 *	Comments:
 *		The invalid memory access in EXIT_UNSUPPORTED/_MSG allow for a stacktrace to be obtained when running the code
 *		using memcheck in valgrind.
 *
 *		EXIT_UNSUPPORTED/EXIT_MSG do not print the "FILE/FUNCTION/LINE" if called before PetscInitilize is called; this
 *		was the original motivation behind the usage of EXIT_BASIC. Now that PetscInitialize is called in main for DTEST
 *		enabled as well, EXIT_BASIC need not be used.
 *
 *	Notation:
 *
 *	References:
 */

#include <stdio.h>
#include <stdlib.h>

#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

#define sign(a) ({ __typeof__ (a) _a = (a); (_a > 0) ? 1 : ((_a < 0) ? -1 : 0); })

#define EXIT_MSG         ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); int *EXIT_VAR = NULL; free(EXIT_VAR); printf("Error: %d\n",EXIT_VAR[1]); exit(1); })
#define EXIT_BASIC       ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); exit(1); })
#define EXIT_UNSUPPORTED ({printf("Error: Unsupported.\n"), EXIT_MSG; })
#define EXIT_ADD_SUPPORT ({printf("Error: Add support.\n"), EXIT_MSG; })
#define EXIT_ERROR(s)    ({printf("Error: %s.\n",s), EXIT_MSG; })
#define PRINT_FILELINE   ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); })
#define FREE_NULL(a)     ({free(a); a = NULL;})
#define UNUSED(x)        (void)(x)


#endif // DPG__Macros_h__INCLUDED
