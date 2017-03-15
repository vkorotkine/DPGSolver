// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Macros_h__INCLUDED
#define DPG__Macros_h__INCLUDED

/*
 *	Purpose:
 *		Define macros.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

#define sign(a) ({ __typeof__ (a) _a = (a); (_a > 0) ? 1 : ((_a < 0) ? -1 : 0); })

#define EXIT_MSG         ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); int *EXIT_VAR = NULL; free(EXIT_VAR); printf("Error: %d\n",EXIT_VAR[1]); exit(1); })
#define EXIT_BASIC       ({ printf("\n\nFILE: %s, FUNCTION: %s (LINE: %d)\n\n\n",__FILE__,__func__,__LINE__); exit(1); })
#define EXIT_UNSUPPORTED ({printf("Error: Unsupported.\n"), EXIT_MSG; })


#endif // DPG__Macros_h__INCLUDED
