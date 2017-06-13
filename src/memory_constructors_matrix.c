// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_constructors_matrix.h"

#include <stdlib.h>

#include "matrix_structs.h"

/*
 *	Purpose:
 *		Provide functions for allocation of matrix structs.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_MATRIX *constructor1_mat (size_t const N0)
{
	struct S_MATRIX *A = calloc(N0 , sizeof *A);
	return A;
}

struct S_MATRIX **constructor2_mat (size_t const N0, size_t const N1)
{
	struct S_MATRIX **A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor1_mat(N1);
	return A;
}


void constructor0_mat_move (char const layout, const char format, const size_t NRows, const size_t NCols,
                            double const *const values, struct S_MATRIX *A)
{
	A = calloc(sizeof *A); // returned

	A->layout = layout;
	A->format = format;
	A->NRows  = NRows;
	A->NCols  = NCols;
	A->values = values;

	return A;
}

void constructor1_mat_move (
