// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "array_free.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "S_OpCSR.h"
#include "matrix_structs.h"

/*
 *	Purpose:
 *		Free dynamically allocated arrays with more than one level of pointer abstraction.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

void free_NULL (void *A)
{
	free(A);
	A = NULL;
}

void array_free2_c(size_t const iMax, char **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_ui(size_t const iMax, unsigned int **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_i(size_t const iMax, int **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_l(size_t const iMax, long **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_ll(size_t const iMax, long long **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_f(size_t const iMax, float **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_d(size_t const iMax, double **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_ld(size_t const iMax, long double **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}

void array_free2_cmplx(size_t const iMax, double complex **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			free_NULL(A[i]);
	free_NULL(A);
}



void array_free3_c(size_t const iMax, size_t const jMax, char ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_c(jMax,A[i]);
	free_NULL(A);
}

void array_free3_ui(size_t const iMax, size_t const jMax, unsigned int ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_ui(jMax,A[i]);
	free_NULL(A);
}

void array_free3_i(size_t const iMax, size_t const jMax, int ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_i(jMax,A[i]);
	free_NULL(A);
}

void array_free3_l(size_t const iMax, size_t const jMax, long ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_l(jMax,A[i]);
	free_NULL(A);
}

void array_free3_ll(size_t const iMax, size_t const jMax, long long ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_ll(jMax,A[i]);
	free_NULL(A);
}

void array_free3_f(size_t const iMax, size_t const jMax, float ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_f(jMax,A[i]);
	free_NULL(A);
}

void array_free3_d(size_t const iMax, size_t const jMax, double ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_d(jMax,A[i]);
	free_NULL(A);
}

void array_free3_ld(size_t const iMax, size_t const jMax, long double ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_ld(jMax,A[i]);
	free_NULL(A);
}

void array_free4_ui(size_t const iMax, size_t const jMax, size_t const kMax, unsigned int ****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free3_ui(jMax,kMax,A[i]);
	free_NULL(A);
}

void array_free4_d(size_t const iMax, size_t const jMax, size_t const kMax, double ****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free3_d(jMax,kMax,A[i]);
	free_NULL(A);
}

void array_free5_d(size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, double *****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free4_d(jMax,kMax,lMax,A[i]);
	free_NULL(A);
}

void array_free1_CSR_d(struct S_OpCSR *A)
{
	if (A) {
		free_NULL(A->rowIndex);
		free_NULL(A->columns);
		free_NULL(A->values);
		free_NULL(A);
	}
}

void array_free4_CSR_d(size_t const iMax, size_t const jMax, size_t const kMax, struct S_OpCSR ****A)
{
	for (size_t i = 0; i < iMax; i++) {
	if (A[i]) {
		for (size_t j = 0; j < jMax; j++) {
		if (A[i][j]) {
			for (size_t k = 0; k < kMax; k++) {
			if (A[i][j][k]) {
				free_NULL(A[i][j][k]->rowIndex);
				free_NULL(A[i][j][k]->columns);
				free_NULL(A[i][j][k]->values);
				free_NULL(A[i][j][k]);
			}}
			free_NULL(A[i][j]);
		}}
		free_NULL(A[i]);
	}}
	free_NULL(A);
}

void array_free5_CSR_d(size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, struct S_OpCSR *****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free4_CSR_d(jMax,kMax,lMax,A[i]);
	free_NULL(A);
}

void matrix_free (struct S_MATRIX *A)
{
	free_NULL(A->values);

	if (A->format == 'S') {
		free_NULL(A->rowIndex);
		free_NULL(A->columns);
	}

	free_NULL(A);
}

void matrix_free2 (size_t const iMax, struct S_MATRIX **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			matrix_free(A[i]);
	free_NULL(A);
}

void matrix_free3 (size_t const N0, size_t const N1, struct S_MATRIX ***A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			matrix_free2(N1,A[i]);
	free_NULL(A);
}

void matrix_free4 (size_t const N0, size_t const N1, size_t const N2, struct S_MATRIX ****A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			matrix_free3(N1,N2,A[i]);
	free_NULL(A);
}

void matrix_free5 (size_t const N0, size_t const N1, size_t const N2, size_t const N3, struct S_MATRIX *****A)
{
	for (size_t i = 0; i < N0; i++)
		if (A[i])
			matrix_free4(N1,N2,N3,A[i]);
	free_NULL(A);
}
