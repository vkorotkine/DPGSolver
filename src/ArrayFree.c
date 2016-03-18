#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

//#include "database.h"
//#include "parameters.h"
//#include "functions.h"

/* 
  Purpose:
    Free dynamically allocated arrays with more than one level of pointer abstraction.

  Comments:

  Notation:

  References:

*/

void ArrayFree2c(int iMax, char **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2i(int iMax, int **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2l(int iMax, long **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2ll(int iMax, long long **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2f(int iMax, float **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2d(int iMax, double **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}

void ArrayFree2ld(int iMax, long double **A) {
  int i;

  for (i = 0; i < iMax; i++) free(A[i]);
  free(A);
}




void ArrayFree3c(int iMax, int jMax, char ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3i(int iMax, int jMax, int ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3l(int iMax, int jMax, long ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3ll(int iMax, int jMax, long long ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3f(int iMax, int jMax, float ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3d(int iMax, int jMax, double ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

void ArrayFree3ld(int iMax, int jMax, long double ***A) {
  int i, j;

  for (i = 0; i < iMax; i++) {
    for (j = 0; j < jMax; j++) free(A[i][j]);
    free(A[i]);
  }
  free(A);
}

