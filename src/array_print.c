// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "array_print.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include "Parameters.h"
#include "Macros.h"

#include "matrix_structs.h"
#include "array_free.h"

/*
 *	Purpose:
 *		Print arrays for debugging purposes.
 *
 *	Comments:
 *		If an elegant way is found to print arbitrary array types, revive array_print (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 *
 */

void array_print_ui(const unsigned int m, const unsigned int n, const unsigned int *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_i(const unsigned int m, const unsigned int n, const int *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 11d ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 11d ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_l(const unsigned int m, const unsigned int n, const long *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12ld ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12ld ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_ll(const unsigned int m, const unsigned int n, const long long *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12lld ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12lld ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_f(const unsigned int m, const unsigned int n, const float *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3f ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3f ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_d(const unsigned int m, const unsigned int n, const double *A, const char layout)
{
	unsigned int i, j;

	bool const AllowZero = 0;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				if ((isnan(A[i*n+j]) || fabs(A[i*n+j]) > EPS) || !AllowZero)
					printf("% .4e ",A[i*n+j]);
				else
					printf(" %d          ",0);
//				printf("% .4e ",A[i*n+j]);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				if ((isnan(A[i+j*m]) || fabs(A[i+j*m]) > EPS) || !AllowZero)
					printf("% .4e ",A[i+j*m]);
				else
					printf(" %d          ",0);
//				printf("% .4e ",A[i+j*m]);
			}
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_ld(const unsigned int m, const unsigned int n, const long double *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3Le ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3Le ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_cmplx(const unsigned int m, const unsigned int n, const double complex *A, const char layout)
{
	unsigned int   i, j;
	double         Ar, Ac;
	double complex Arc;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				Arc = A[i*n+j];
				Ar   = creal(Arc);
				Ac   = cimag(Arc);
				printf("% .4e%c%.4ei ",Ar,(Ac >= 0.0)?'+':'\0',Ac);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				Arc = A[i+j*m];
				Ar   = creal(Arc);
				Ac   = cimag(Arc);
				printf("% .4e%c%.4ei ",Ar,(Ac >= 0.0)?'+':'\0',Ac);
			}
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void matrix_print (struct S_MATRIX const *const A)
{
	array_print_d(A->extents[0],A->extents[1],A->data,A->layout);
}

void multiarray_print (struct S_MULTI_ARRAY const *const A)
{
	size_t const *const extents = A->extents;

	printf("Multi-array extents:");
	for (size_t i = 0; i < A->order; i++)
		printf(" %zu",A->extents[i]);
	printf("\n\n");

	switch (A->order) {
	case 2: {
		struct S_MATRIX const A_M = constructor_matrix1_move_multiarray2_2(A);
		printf("(:,:):\n");
		matrix_print(&A_M);
		break;
	} case 3:
		for (size_t i = 0; i < extents[2]; i++) {
			struct S_MATRIX const A_M = constructor_matrix1_move_multiarray3_2(A,i);

			printf("(:,:,%zu):\n",i);
			matrix_print(&A_M);
		}
		break;
	case 4:
		for (size_t i = 0; i < extents[2]; i++) {
		for (size_t j = 0; j < extents[3]; j++) {
			struct S_MATRIX const A_M = constructor_matrix1_move_multiarray4_2(A,i,j);

			printf("(:,:,%zu,%zu):\n",i,j);
			matrix_print(&A_M);
		}}
		break;
	default:
		printf("%zu\n",A->order);
		EXIT_UNSUPPORTED;
		break;
	}
}
