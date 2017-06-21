// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "matrix_structs.h"

#include <stdlib.h>

#include "array_free.h"

/*
 *	Purpose:
 *		Provide functions relating to matrix structs.
 *
 *	Comments:
 *		The 'pointer' constructors can only return pointers to pointers to structs and never the actual pointers to the
 *		structs at the end of the chain. These functions are used to initialize the storage for the potential operators
 *		which are then computed in setup_operators_*.
 *
 *	Notation:
 *		The function naming convention for matrix constructors is:
 *
 *		constructor_[0][1]_[2_3]_[4][5]:
 *			[0]: Type of struct to be constructed.
 *				Options : matrix, multiarray
 *			[1]: Number indicating the level of dereferencing in the returned matrix struct.
 *			[2_3]: Naming convention for type of constructor and from what kind of data.
 *				Options [2]           : default, pointer, empty, copy, move
 *				        [3] (Optional): mat (default if not specified), (d)ouble, multiarray
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


// Pointer constructors (minimum dereferencing of 2)

struct S_MATRIX **constructor_matrix2_pointer (size_t const N0)
{
	struct S_MATRIX **A = calloc(N0 , sizeof *A);
	return A;
}

struct S_MATRIX ***constructor_matrix3_pointer (size_t const N0, size_t const N1)
{
	struct S_MATRIX ***A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor_matrix2_pointer(N1);
	return A;
}

struct S_MATRIX ****constructor_matrix4_pointer (size_t const N0, size_t const N1, size_t const N2)
{
	struct S_MATRIX ****A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor_matrix3_pointer(N1,N2);
	return A;
}

struct S_MATRIX *****constructor_matrix5_pointer (size_t const N0, size_t const N1, size_t const N2, size_t const N3)
{
	struct S_MATRIX *****A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor_matrix4_pointer(N1,N2,N3);
	return A;
}


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

struct S_MATRIX *constructor_matrix1_empty (char const layout, size_t const m, size_t const n)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a (D)ense matrix struct where data memory has been allocated.
	 */

	struct S_MATRIX *A = constructor_matrix1_default(); // returned

	A->layout     = layout;
	A->extents[0] = m;
	A->extents[1] = n;
	A->data       = calloc(m*n , sizeof *(A->data)); // keep

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

struct S_MULTI_ARRAY *constructor_multiarray1_empty3 (char const layout, size_t const N0, size_t const N1,
                                                      size_t const N2)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a 3D multi-array struct where data memory has been allocated.
	 */

	size_t const order = 3;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = layout;
	A->extents[0] = N0;
	A->extents[1] = N1;
	A->extents[2] = N2;
	A->data       = calloc(N0*N1*N2 , sizeof *(A->data)); // keep

	return A;
}

struct S_MULTI_ARRAY *constructor_multiarray1_empty4 (char const layout, size_t const N0, size_t const N1,
                                                      size_t const N2, size_t const N3)
{
	/*
	 *	Purpose:
	 *		Return a pointer to a 4D multi-array struct where data memory has been allocated.
	 */

	size_t const order = 4;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = layout;
	A->extents[0] = N0;
	A->extents[1] = N1;
	A->extents[2] = N2;
	A->extents[3] = N3;
	A->data       = calloc(N0*N1*N2*N3 , sizeof *(A->data)); // keep

	return A;
}


// Copy constructors

struct S_MATRIX *constructor_matrix1_copy (struct S_MATRIX const *const A)
{
	struct S_MATRIX *B = constructor_matrix1_default(); // returned

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

struct S_MATRIX constructor_matrix1_move_multiarray2_2 (struct S_MULTI_ARRAY const *const A)
{
	struct S_MATRIX B;

	size_t const order = 2;
	for (size_t i = 0; i < order; i++)
		B.extents[i] = A->extents[i];

	B.layout = A->layout;
	B.data   = A->data;

	return B;
}

struct S_MATRIX constructor_matrix1_move_multiarray3_2 (struct S_MULTI_ARRAY const *const A, size_t const N2)
{

	struct S_MATRIX B;

	size_t const order = 2;
	for (size_t i = 0; i < order; i++)
		B.extents[i] = A->extents[i];

	B.layout = A->layout;

	size_t const size = compute_size(order,B.extents);
	B.data = &A->data[size*(N2)];

	return B;
}

struct S_MATRIX constructor_matrix1_move_multiarray4_2 (struct S_MULTI_ARRAY const *const A, size_t const N2,
                                                        size_t const N3)
{
	struct S_MATRIX B;

	size_t const order = 2;
	for (size_t i = 0; i < order; i++)
		B.extents[i] = A->extents[i];

	B.layout = A->layout;

	size_t const size = compute_size(order,B.extents);
	B.data = &A->data[size*(N2+(A->extents[2])*(N3))];

	return B;
}

struct S_MATRIX constructor_matrix1_move_multiarray4_1 (struct S_MULTI_ARRAY const *const A, size_t const N1,
                                                        size_t const N2, size_t const N3)
{
	struct S_MATRIX B;
	B.structure = oneD_M;

	size_t const order = 1;
	for (size_t i = 0; i < order; i++)
		B.extents[i] = A->extents[i];

	B.layout = A->layout;

	size_t const size = compute_size(order,B.extents);
	B.data = &A->data[size*(N1+(A->extents[1])*(N2+(A->extents[2])*(N3)))];

	return B;
}

struct S_MULTI_ARRAY *constructor_multiarray1_move_d_2 (char const layout, size_t const N0, size_t const N1,
                                                        double *const data)
{
	/*
	 *	Purpose:
	 *		Move (d)ouble to multiarray of order 2.
	 */

	size_t const order = 2;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = layout;
	A->extents[0] = N0;
	A->extents[1] = N1;
	A->data       = data;

	return A;
}

struct S_MULTI_ARRAY *constructor_multiarray1_move_d_3 (char const layout, size_t const N0, size_t const N1,
                                                        size_t const N2, double *const data)
{
	/*
	 *	Purpose:
	 *		Move (d)ouble to multiarray of order 3.
	 */

	size_t const order = 3;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = layout;
	A->extents[0] = N0;
	A->extents[1] = N1;
	A->extents[2] = N2;
	A->data       = data;

	return A;
}

struct S_MULTI_ARRAY *constructor_multiarray1_move_d_4 (char const layout, size_t const N0, size_t const N1,
                                                        size_t const N2, size_t const N3, double *const data)
{
	/*
	 *	Purpose:
	 *		Move (d)ouble to multiarray of order 4.
	 */

	size_t const order = 4;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = layout;
	A->extents[0] = N0;
	A->extents[1] = N1;
	A->extents[2] = N2;
	A->extents[3] = N3;
	A->data       = data;

	return A;
}

struct S_MULTI_ARRAY *constructor_multiarray1_move_matrix_2 (struct S_MATRIX const *const A_M, bool const destruct_M)
{
	/*
	 *	Purpose:
	 *		Move information from a matrix to a multiarray struct.
	 *
	 *	Comments:
	 *		the destruct_M parameter indicates whether the "shell" of the passed matrix struct should be freed.
	 */
	size_t const order = 2;
	struct S_MULTI_ARRAY *A = constructor_multiarray1_default(order); // returned

	A->layout     = A_M->layout;
	A->extents[0] = A_M->extents[0];
	A->extents[1] = A_M->extents[1];
	A->data       = A_M->data;

	if (destruct_M)
		destructor_matrix1_default_const(A_M);

	return A;
}


// destructors

void destructor_matrix1_default (struct S_MATRIX *A)
{
	/*
	 *	Comments:
	 *		Does not free A->data.
	 */

	free_NULL(A);
}

void destructor_matrix1_default_const (struct S_MATRIX const *const A)
{
	destructor_matrix1_default((void *) A);
}

void destructor_matrix1 (struct S_MATRIX *A)
{
	free_NULL(A->data);
	free_NULL(A);
}

void destructor_matrix2_pointer (size_t const N0, struct S_MATRIX **A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			destructor_matrix1(A[i]);
	free_NULL(A);
}

void destructor_matrix3_pointer (size_t const N0, size_t const N1, struct S_MATRIX ***A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			destructor_matrix2_pointer(N1,A[i]);
	free_NULL(A);
}

void destructor_matrix4_pointer (size_t const N0, size_t const N1, size_t const N2, struct S_MATRIX ****A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			destructor_matrix3_pointer(N1,N2,A[i]);
	free_NULL(A);
}

void destructor_matrix5_pointer (size_t const N0, size_t const N1, size_t const N2, size_t const N3,
                                 struct S_MATRIX *****A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			destructor_matrix4_pointer(N1,N2,N3,A[i]);
	free_NULL(A);
}


void destructor_multiarray1_default (struct S_MULTI_ARRAY *A)
{
	free_NULL(A->extents);
	free_NULL(A);
}

void destructor_multiarray1_default_const (struct S_MULTI_ARRAY const *const A)
{
	destructor_multiarray1_default((void *) A);
}

void destructor_multiarray1 (struct S_MULTI_ARRAY *A)
{
	free_NULL(A->extents);
	free_NULL(A->data);
	free_NULL(A);
}

void destructor_multiarray1_const (struct S_MULTI_ARRAY const *const A)
{
	destructor_multiarray1((void *) A);
}


// Additional functions

static inline size_t compute_size (size_t const NExt, size_t const *const extents)
{
	size_t size = 1;
	for (size_t i = 0; i < NExt; i++)
		size *= extents[i];

	return size;
}

void set_to_zero_multiarray (struct S_MULTI_ARRAY *const A)
{
	double *A_data = A->data;
	for (size_t size = compute_size(A->order,A->extents); size--; )
		*A_data++ = 0.0;
}
