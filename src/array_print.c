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

void array_print_i(int m, int n, int *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12d ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void array_print_l(int m, int n, long *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12ld ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void array_print_ll(int m, int n, long long *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% 12lld ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void array_print_f(int m, int n, float *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3f ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void array_print_d(int m, int n, double *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3e ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void array_print_ld(int m, int n, long double *A) {
  int i, j;

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf("% .3Le ",A[i*n+j]);
    }
    printf("\n");
  }
  printf("\n");
}
