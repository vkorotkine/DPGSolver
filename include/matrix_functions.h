// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__matrix_functions_h__INCLUDED
#define DPG__matrix_functions_h__INCLUDED

#include "mkl.h"

#include "S_OpCSR.h"
#include "matrix_structs.h"


extern double *diag_d          (const double *x, const unsigned int N);
extern double *identity_d      (const unsigned int N);
extern double *inverse_d       (const unsigned int N, const unsigned int NRHS, const double *A, const double *b);
extern void   mm_diag_d        (const unsigned int NRows, const unsigned int NCols, double const *const a,
                                double const *const A, double *const Output, const double alpha, const double beta,
                                const char side, const char layout);
extern double *mm_Alloc_d      (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double *A, const double *B);
extern void   mm_d             (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double beta,
                                const double *A, const double *B, double *C);
extern void   mm_dcc           (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                                const int m, const int n, const int k, const double alpha, const double beta,
                                const double *A, const void *B, void *C);
extern void   mm_CTN_d         (const int m, const int n, const int k, const double *A, const double *B, double *C);
extern void   mm_CTN_CSR_d     (const int m, const int n, const int k, const struct S_OpCSR *A, const double *B, double *C);
extern void   convert_to_CSR_d (const unsigned int NRows, const unsigned int NCols, const double *Input,
                                struct S_OpCSR **Output);

extern struct S_MATRIX *mm_matrix_alloc (char const layout, char const opA, char const opB, struct S_MATRIX const *const A,
                                         struct S_MATRIX const *const B);
extern void mm_matrix (char const opA, char const opB, double const alpha, double const beta,
                       struct S_MATRIX const *const A, struct S_MATRIX const *const B, struct S_MATRIX *const C);
extern struct S_MATRIX *inverse_mat  (struct S_MATRIX const *const A);
extern struct S_MATRIX *identity_mat (unsigned int const N);
extern struct S_MATRIX *mat_copy     (struct S_MATRIX const *const A);
extern struct S_MATRIX *mm_diag_matrix_alloc (char const layoutB, char const opA, char const side,
                                              struct S_MATRIX const *const a, struct S_MATRIX const *const A);
extern void mm_diag_matrix (char const opA, char const side, const double alpha, const double beta,
                            struct S_MATRIX const *const a, struct S_MATRIX const *const A, struct S_MATRIX *const B);
extern void transpose_matrix   (struct S_MATRIX *const A);
extern void set_to_zero_matrix (struct S_MATRIX *const A);


extern struct S_MATRIX *mat_constructor_move (char const layout, char const format, size_t const NRows,
                                              size_t const NCols, double *const values, struct S_OpCSR *const A_CSR);
extern struct S_MATRIX **mat_constructor2_move (char const layout, char const format, size_t const NRows,
                                                size_t const NCols, size_t const dim1, double *const *const values,
                                                struct S_OpCSR *const *const A_CSR);


#endif // DPG__matrix_functions_h__INCLUDED
