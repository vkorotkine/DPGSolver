// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__matrix_functions_h__INCLUDED
#define DPG__matrix_functions_h__INCLUDED

#include "mkl.h"

#include "S_OpCSR.h"


extern double *diag_d          (const double *x, const unsigned int N);
extern double *identity_d      (const unsigned int N);
extern double *inverse_d       (const unsigned int N, const unsigned int NRHS, const double *A, const double *b);
extern void   mm_diag_d        (const unsigned int NRows, const unsigned int NCols, double *a, double *A, double *Output,
                                const double alpha, const char side, const char layout);
extern double *mm_Alloc_d      (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double *A, const double *B);
extern void   mm_d             (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double beta,
                                const double *A, const double *B, double *C);
extern void   mm_dcc           (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double beta,
                                double *A, void *B, void *C);
extern void   mm_CTN_d         (const int m, const int n, const int k, double *A, double *B, double *C);
extern void   mm_CTN_CSR_d     (const int m, const int n, const int k, const struct S_OpCSR *A, double *B, double *C);
extern void   convert_to_CSR_d (const unsigned int NRows, const unsigned int NCols, const double *Input,
                                struct S_OpCSR **Output);

#endif // DPG__matrix_functions_h__INCLUDED
