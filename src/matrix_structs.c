// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "matrix_structs.h"

/*
 *	Purpose:
 *		Provide functions relating to matrix structs.
 *
 *	Notation:
 *		The function naming convention for matrix constructors is:
 *
 *		constructor_mat[0]_[1_2]_[3][4]:
 *			[0]: Number indicating the level of dereferencing in the returned matrix struct.
 *			[1_2]: Naming convention for type of constructor and from what kind of data.
 *				Options [1]           : default, empty, copy, move
 *				        [2] (Optional): mat (default if not specified), (d)ouble
 *			[3]: Matrix format parameter.
 *				Options: 'D'ense, 'S'parse.
 *			[4]: Matrix order parameter (i.e. the dimension of the matrix)
 */

static inline size_t compute_size (size_t const NExt, size_t const *const extent);


// Empty constructors

struct S_MATRIX *constructor_mat1_default (const size_t order)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a matrix struct where only the extent (and not the data) has been allocated.
	 */

	struct S_MATRIX *A = calloc(1 , sizeof *A); // returned;

	A->extent = calloc(order , sizeof *(A->extent)); // keep
	A->order  = order;

	if (order == 1) {
		// Vectors are dense columns
		A->layout = 'C';
		A->format = 'D';
	}

	return A;
}

struct S_MATRIX *constructor_mat1_empty_D2 (const char layout, const size_t N0, const size_t N1)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a (D)ense (2)-dimensional matrix struct where both the extent and the data memory have
	 *		been allocated.
	 *
	 *	Notation:
	 *		layout : 'R'ow or 'C'olumn major
	 *		N0     : (N)umber of Rows
	 *		N1     : (N)umber of Cols
	 */

	struct S_MATRIX *A = calloc(1 , sizeof *A); // returned

	A->format = 'D';
	A->extent = malloc(2 * sizeof *(A->extent)); // keep

	A->layout = layout;
	A->extent[0] = N0;
	A->extent[1] = N1;
	A->values = calloc(N0*N1 , sizeof *(A->values)); // keep

	return A;
}


// Copy constructors

struct S_MATRIX *constructor_mat1_copy_D (struct S_MATRIX const *const A)
{
	struct S_MATRIX *B = calloc(1 , sizeof *B); // returned

	size_t const order = A->order;
	B->extents = malloc(order * sizeof *(B->extents)); // keep

	for (size_t i = 0; i < order; i++)
		B->extents[i] = A->extents[i];

	B->layout = A->layout;
	B->format = A->format;
	B->order  = A->order;

	size_t const size = compute_size(order,B->extents);

	B->data = malloc(size * sizeof *(B->data)); // keep
	for (size_t i = 0; i < size; i++)
		B->data[i] = A->data[i];

	return B;
}


// Move constructors

struct S_MATRIX *constructor_mat1_move_d_D1 (size_t const N0, double *const data)
{
	size_t const order = 1;

	struct S_MATRIX *a = constructor_mat1_default(order);

	a->extent[0] = N0;
	a->data      = data;
}




// destructors





// Additional functions

static inline size_t compute_size (size_t const NExt, size_t const *const extent)
{
	size_t size = 1;
	for (size_t i = 0; i < NExt; i++)
		size *= extent[i];

	return size;
}
