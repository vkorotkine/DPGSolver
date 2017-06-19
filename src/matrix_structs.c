// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "matrix_structs.h"

#include <stdlib.h>

/*
 *	Purpose:
 *		Provide functions relating to matrix structs.
 *
 *	Notation:
 *		The function naming convention for matrix constructors is:
 *
 *		constructor_[0][1]_[2_3]_[4][5]:
 *			[0]: Type of struct to be constructed.
 *				Options : matrix, multiarray
 *			[1]: Number indicating the level of dereferencing in the returned matrix struct.
 *			[2_3]: Naming convention for type of constructor and from what kind of data.
 *				Options [2]           : default, empty, copy, move
 *				        [3] (Optional): mat (default if not specified), (d)ouble
 *			[4] (Optional): Matrix format parameter.
 *				Options: 'D'ense (default), 'S'parse.
 *			[5]: 'order' of the resulting struct.
 *				Options : 2 (default), 1 (matrix (vector)), other
 *
 *		For the move_multiarray functions, the notation is similar to above with the two trailing numbers indicating the
 *		dimension of the input multi_array and the output matrix. For example 'constructor_matrix1_move_multiarray41'
 *		returns a S_MATRIX* of order 1 from an input S_MULTI_ARRAY* of order 4 where the data array points to the
 *		appropriate location in memory based on multi-dimensional array indexing.
 */

static inline size_t compute_size (size_t const NExt, size_t const *const extents);


// Empty constructors

struct S_MATRIX *constructor_matrix1_default (void)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a matrix struct.
	 */

	struct S_MATRIX *A = calloc(1 , sizeof *A); // returned
	return A;
}

struct S_MULTI_ARRAY *constructor_multiarray1_default (const size_t order)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a multi-array struct where only the extents (and not the data) has been allocated.
	 */

	struct S_MULTI_ARRAY *A = calloc(1 , sizeof *A); // returned

	A->extents = calloc(order , sizeof *(A->extents)); // keep
	A->order   = order;

	return A;
}

struct S_MATRIX *constructor_matrix1_empty (char const layout, size_t const m, size_t const n)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a (D)ense matrix struct where data memory has been allocated.
	 */

	struct S_MATRIX *A = calloc(1 , sizeof *A); // returned

	A->layout     = layout;
	A->extents[0] = m;
	A->extents[1] = n;
	A->data       = calloc(m*n , sizeof *(A->data)); // keep

	return A;
}


// Copy constructors

struct S_MATRIX *constructor_matrix1_copy (struct S_MATRIX const *const A)
{
	struct S_MATRIX *B = calloc(1 , sizeof *B); // returned

	for (size_t i = 0; i < 2; i++)
		B->extents[i] = A->extents[i];

	B->layout = A->layout;

	size_t const size = compute_size(2,B->extents);

	B->data = malloc(size * sizeof *(B->data)); // keep
	for (size_t i = 0; i < size; i++)
		B->data[i] = A->data[i];

	return B;
}


// Move constructors

struct S_MATRIX *constructor_matrix1_move_d_1 (size_t const m, double *const data)
{
	struct S_MATRIX *a = constructor_matrix1_default(); // returned

	a->structure  = oneD_M;
	a->extents[0] = m;
	a->data       = data;

	return a;
}



// Move constructors from multi-arrays

struct S_MATRIX *constructor_matrix1_move_multiarray32(struct S_MULTI_ARRAY *A, size_t const N2)
{
	size_t const order = 2;

	struct S_MATRIX *B = constructor_matrix1_default(); // returned

	for (size_t i = 0; i < order; i++)
		B->extents[i] = A->extents[i];

	B->layout = A->layout;

	size_t const size = compute_size(order,B->extents);
	B->data = &A->data[size*(N2)];

	return B;
}

struct S_MATRIX *constructor_matrix1_move_multiarray42(struct S_MULTI_ARRAY *A, size_t const N2, size_t const N3)
{
	size_t const order = 2;

	struct S_MATRIX *B = constructor_matrix1_default(); // returned

	for (size_t i = 0; i < order; i++)
		B->extents[i] = A->extents[i];

	B->layout = A->layout;

	size_t const size = compute_size(order,B->extents);
	B->data = &A->data[size*(N2+(A->extents[2])*(N3))];

	return B;
}

struct S_MATRIX *constructor_matrix1_move_multiarray41(struct S_MULTI_ARRAY *A, size_t const N2, size_t const N3,
                                                       size_t const N4)
{
	size_t const order = 1;

	struct S_MATRIX *B = constructor_matrix1_default(); // returned
	B->structure = oneD_M;

	for (size_t i = 0; i < order; i++)
		B->extents[i] = A->extents[i];

	B->layout = A->layout;

	size_t const size = compute_size(order,B->extents);
	B->data = &A->data[size*(N2+(A->extents[2])*(N3+(A->extents[3])*(N4)))];

	return B;
}



// destructors





// Additional functions

static inline size_t compute_size (size_t const NExt, size_t const *const extents)
{
	size_t size = 1;
	for (size_t i = 0; i < NExt; i++)
		size *= extents[i];

	return size;
}
