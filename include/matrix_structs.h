// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_structs_h__INCLUDED
#define DPG__matrix_structs_h__INCLUDED

#include <stddef.h>

struct S_MATRIX {
	char   layout, // May be (R)ow or (C)olumn major.
	       format; // May be (D)ense or (S)parse.
	size_t NRows, NCols;

	double *values;

	// Used for CSR format
	unsigned int *rowIndex, *columns;
};

#endif // DPG__matrix_structs_h__INCLUDED
