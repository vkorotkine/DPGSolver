// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_structs_h__INCLUDED
#define DPG__matrix_structs_h__INCLUDED

#include <stddef.h>

struct S_MATRIX {
	char   layout,  // May be (R)ow or (C)olumn major.
	       format;  // May be (D)ense or (S)parse.

	size_t order,   // Number of matrix dimensions (== 2 for std matrices)
	       *extent; // Size of arrays in each dimension
//	       size;    // Total number of entries in the matrix

	double *data;

	// Used for CSR format
	unsigned int *rowIndex, *columns;
};

extern struct S_MATRIX *constructor1_default_mat (const size_t order);

#endif // DPG__matrix_structs_h__INCLUDED
