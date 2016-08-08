// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__macros_h__INCLUDED
#define DPG__macros_h__INCLUDED

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

#define EXIT_MSG ({ printf("FILE: %s, FUNCTION: %s (LINE: %d)\n",__FILE__,__func__,__LINE__); exit(1); })


#endif // DPG__macros_h__INCLUDED
