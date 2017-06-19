// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_structs_h__INCLUDED
#define DPG__matrix_structs_h__INCLUDED

#include <stdbool.h>
#include <stddef.h>

struct S_MATRIX {
	/*
	 *	Purpose:
	 *		Defines a struct to support dense (double) matrices (2D).
	 *
	 *	Comments:
	 *		Primarily used to store operators.
	 *		Can be used as a 1D vector by defining structure = oneD_M.
	 */

	enum matrix_structure { full_M, diag_M, identity_M, oneD_M } structure;

	char layout; // May be (R)ow or (C)olumn major.

	size_t extents[2]; // Number of rows/columns

	double *data;
};

extern struct S_MATRIX *constructor_matrix1_default (void);
extern struct S_MATRIX *constructor_matrix1_copy (struct S_MATRIX const *const A);

struct S_MULTI_ARRAY {
	/*
	 *	Purpose:
	 *		Defines a struct to support dense multi-dimensional (double) arrays.
	 *
	 *	Comments:
	 *		The multi-array struct is intended to be used as a higher-dimensional matrix where move constructors are
	 *		used to form matrix structs for appropriate sub-blocks. As the data is stored contiguously in memory, the
	 *		multi-array may also be acted on over multiple dimensions at once.
	 */

	char layout; // May be (R)ow or (C)olumn major.

	size_t order,    // Number of dimensions.
	       *extents; // Size of arrays in each dimension

	double *data;
};

#endif // DPG__matrix_structs_h__INCLUDED
