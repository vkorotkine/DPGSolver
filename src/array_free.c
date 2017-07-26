// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "array_free.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "Macros.h"
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

void array_free2_c(size_t const iMax, char **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			FREE_NULL(A[i]);
	FREE_NULL(A);
}

void array_free2_ui(size_t const iMax, unsigned int **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			FREE_NULL(A[i]);
	FREE_NULL(A);
}

void array_free2_i(size_t const iMax, int **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			FREE_NULL(A[i]);
	FREE_NULL(A);
}

void array_free2_d(size_t const iMax, double **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			FREE_NULL(A[i]);
	FREE_NULL(A);
}

void array_free2_cmplx(size_t const iMax, double complex **A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			FREE_NULL(A[i]);
	FREE_NULL(A);
}



void array_free3_c(size_t const iMax, size_t const jMax, char ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_c(jMax,A[i]);
	FREE_NULL(A);
}

void array_free3_ui(size_t const iMax, size_t const jMax, unsigned int ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_ui(jMax,A[i]);
	FREE_NULL(A);
}

void array_free3_i(size_t const iMax, size_t const jMax, int ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_i(jMax,A[i]);
	FREE_NULL(A);
}

void array_free3_d(size_t const iMax, size_t const jMax, double ***A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free2_d(jMax,A[i]);
	FREE_NULL(A);
}



void array_free4_ui(size_t const iMax, size_t const jMax, size_t const kMax, unsigned int ****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free3_ui(jMax,kMax,A[i]);
	FREE_NULL(A);
}

void array_free4_d(size_t const iMax, size_t const jMax, size_t const kMax, double ****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free3_d(jMax,kMax,A[i]);
	FREE_NULL(A);
}

void array_free5_d(size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, double *****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free4_d(jMax,kMax,lMax,A[i]);
	FREE_NULL(A);
}

void array_free1_CSR_d(struct S_OpCSR *A)
{
	if (A) {
		FREE_NULL(A->rowIndex);
		FREE_NULL(A->columns);
		FREE_NULL(A->values);
		FREE_NULL(A);
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
				FREE_NULL(A[i][j][k]->rowIndex);
				FREE_NULL(A[i][j][k]->columns);
				FREE_NULL(A[i][j][k]->values);
				FREE_NULL(A[i][j][k]);
			}}
			FREE_NULL(A[i][j]);
		}}
		FREE_NULL(A[i]);
	}}
	FREE_NULL(A);
}

void array_free5_CSR_d(size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, struct S_OpCSR *****A)
{
	for (size_t i = 0; i < iMax; i++)
		if (A[i])
			array_free4_CSR_d(jMax,kMax,lMax,A[i]);
	FREE_NULL(A);
}
