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

void array_free2_c(int iMax, char **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_i(int iMax, int **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_l(int iMax, long **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_ll(int iMax, long long **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_f(int iMax, float **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_d(int iMax, double **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}

void array_free2_ld(int iMax, long double **A)
{
	int i;

	for (i = 0; i < iMax; i++)
		free(A[i]);
	free(A);
}



void array_free3_c(int iMax, int jMax, char ***A)
{
	int i, j;

	for (i = 0; i < iMax; i++) {
		for (j = 0; j < jMax; j++)
			free(A[i][j]);
		free(A[i]);
	}
	free(A);
}

void array_free3_i(int iMax, int jMax, int ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void array_free3_l(int iMax, int jMax, long ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void array_free3_ll(int iMax, int jMax, long long ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void array_free3_f(int iMax, int jMax, float ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void array_free3_d(int iMax, int jMax, double ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void array_free3_ld(int iMax, int jMax, long double ***A)
{
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++)
		free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

