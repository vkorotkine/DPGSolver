// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "array_print.h"

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

void array_print_i(const unsigned int m, const unsigned int n, int *A, const char layout)
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

void array_print_l(const unsigned int m, const unsigned int n, long *A, const char layout)
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

void array_print_ll(const unsigned int m, const unsigned int n, long long *A, const char layout)
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

void array_print_f(const unsigned int m, const unsigned int n, float *A, const char layout)
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

void array_print_d(const unsigned int m, const unsigned int n, double *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .4e ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .4e ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_ld(const unsigned int m, const unsigned int n, long double *A, const char layout)
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
