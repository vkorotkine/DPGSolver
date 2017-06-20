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


extern struct S_MATRIX **   constructor_matrix2_pointer (size_t const N0);
extern struct S_MATRIX ***  constructor_matrix3_pointer (size_t const N0, size_t const N1);
extern struct S_MATRIX **** constructor_matrix4_pointer (size_t const N0, size_t const N1, size_t const N2);
extern struct S_MATRIX *****constructor_matrix5_pointer (size_t const N0, size_t const N1, size_t const N2,
                                                         size_t const N3);

extern struct S_MATRIX *constructor_matrix1_default (void);
extern struct S_MATRIX *constructor_matrix1_copy (struct S_MATRIX const *const A);

extern struct S_MATRIX *constructor_matrix1_move_d_1 (size_t const m, double *const data);
extern struct S_MATRIX constructor_matrix1_move_multiarray2_2 (struct S_MULTI_ARRAY const *const A);
extern struct S_MATRIX constructor_matrix1_move_multiarray3_2 (struct S_MULTI_ARRAY const *const A, size_t const N2);
extern struct S_MATRIX constructor_matrix1_move_multiarray4_2 (struct S_MULTI_ARRAY const *const A, size_t const N2,
                                                               size_t const N3);
extern struct S_MATRIX constructor_matrix1_move_multiarray4_1 (struct S_MULTI_ARRAY const *const A, size_t const N1,
	                                                           size_t const N2, size_t const N3);
extern void destructor_matrix1_default (struct S_MATRIX *A);
extern void destructor_matrix1_default_const (struct S_MATRIX const *const A);


extern struct S_MULTI_ARRAY *constructor_multiarray1_empty4 (char const layout, size_t const N0, size_t const N1,
                                                             size_t const N2, size_t const N3);

extern struct S_MULTI_ARRAY *constructor_multiarray1_move_d_2 (char const layout, size_t const N0, size_t const N1,
                                                               double *const data);
extern struct S_MULTI_ARRAY *constructor_multiarray1_move_d_3 (char const layout, size_t const N0, size_t const N1,
	                                                           size_t const N2, double *const data);
extern struct S_MULTI_ARRAY *constructor_multiarray1_move_d_4 (char const layout, size_t const N0, size_t const N1,
                                                               size_t const N2, size_t const N3, double *const data);
extern struct S_MULTI_ARRAY *constructor_multiarray1_move_matrix_2 (struct S_MATRIX const *const A_M);
extern void destructor_multiarray1_default (struct S_MULTI_ARRAY *A);
extern void destructor_multiarray1_default_const (struct S_MULTI_ARRAY const *const A);


extern void set_to_zero_multiarray (struct S_MULTI_ARRAY *const A);

#endif // DPG__matrix_structs_h__INCLUDED
