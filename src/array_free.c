#include <stdlib.h>

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

void array_free2_c(unsigned int iMax, char **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ui(unsigned int iMax, unsigned int **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_i(unsigned int iMax, int **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_l(unsigned int iMax, long **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ll(unsigned int iMax, long long **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_f(unsigned int iMax, float **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_d(unsigned int iMax, double **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ld(unsigned int iMax, long double **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}



void array_free3_c(unsigned int iMax, unsigned int jMax, char ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ui(unsigned int iMax, unsigned int jMax, unsigned int ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_i(unsigned int iMax, unsigned int jMax, int ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_l(unsigned int iMax, unsigned int jMax, long ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ll(unsigned int iMax, unsigned int jMax, long long ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_f(unsigned int iMax, unsigned int jMax, float ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_d(unsigned int iMax, unsigned int jMax, double ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ld(unsigned int iMax, unsigned int jMax, long double ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

