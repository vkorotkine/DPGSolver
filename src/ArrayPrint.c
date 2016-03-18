#include <stdlib.h>
#include <stdio.h>
//#include <string.h>

//#include "database.h"
//#include "parameters.h"
//#include "functions.h"

/* 
  Purpose:
    Print arrays for debugging.

  Comments:

  Notation:

  References:

*/

void ArrayPrinti(int m, int n, int *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12d ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void ArrayPrintl(int m, int n, long *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12ld ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void ArrayPrintll(int m, int n, long long *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12lld ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void ArrayPrintf(int m, int n, float *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3f ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void ArrayPrintd(int m, int n, double *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3e ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void ArrayPrintld(int m, int n, long double *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3Le ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}
