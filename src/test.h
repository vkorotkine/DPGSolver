// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_h__INCLUDED
#define DPG__test_h__INCLUDED

/*
 *	Purpose:
 *		Define global functions and structures for testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

// Support functions
extern void test_print (const unsigned int pass);

// Implementation tests
extern void test_imp_array_find_index          (void);
extern void test_imp_array_norm                (void);
extern void test_imp_array_sort                (void);
extern void test_imp_array_swap                (void);
extern void test_imp_math_factorial            (void);
extern void test_imp_math_gamma                (void);
extern void test_imp_matrix_diag               (void);
extern void test_imp_matrix_identity           (void);
extern void test_imp_matrix_inverse            (void);
extern void test_imp_matrix_mm                 (void);
extern void test_imp_find_periodic_connections (void);
extern void test_imp_cubature_TP               (void);
extern void test_imp_cubature_SI               (void);
extern void test_imp_cubature_PYR              (void);
extern void test_imp_basis_TP                  (void);
extern void test_imp_basis_SI                  (void);
extern void test_imp_basis_PYR                 (void);
extern void     get_scaling_basis_TRI          (const unsigned int i, const unsigned int j, const double b,
                                                double *con_i, double *con_j, double *con_b);
extern void     get_scaling_basis_TET          (const unsigned int i, const unsigned int j, const unsigned int k,
                                                const double b, const double c, double *con_i, double *con_j,
                                                double *con_k, double *con_b, double *con_c);
extern void     get_scaling_basis_PYR          (const unsigned int i, const unsigned int j, const unsigned int k,
                                                const double c, double *con_i, double *con_j, double *con_k,
                                                double *con_c);
extern void     poly2                          (const double *r, const double *s, const double *t,
                                                const unsigned int Nn, double **f, double **f_r, double **f_s,
                                                double **f_t);
extern void test_imp_grad_basis_TP             (void);
extern void test_imp_grad_basis_SI             (void);
extern void test_imp_grad_basis_PYR            (void);
extern void test_imp_plotting                  (void);
extern void test_imp_sum_factorization         (void);
extern void test_imp_fluxes_inviscid           (void);
extern void test_imp_get_facet_ordering        (void);

// Speed tests
extern void test_speed_array_swap (void);
extern void     array_swap1_ui    (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
                                   const unsigned int stepIn);
extern void     array_swap2_ui    (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
                                   const unsigned int stepIn);
extern void     array_swap1_d     (register double *arr1, register double *arr2, const unsigned int NIn,
                                   const unsigned int stepIn);
extern void     array_swap2_d     (register double *arr1, register double *arr2, const unsigned int NIn,
                                   const unsigned int stepIn);
extern void test_speed_mm_CTN     (void);
extern void     mm_CTN_mv_fully_unrolled_d (const int m, const int n, const int k, double *A, double *B, double *C);
extern void     mm_CTN_mvBLAS_d            (const int m, const int n, const int k, double *A, double *B, double *C);
extern void     mm_CTN_fully_unrolled_mv_d (const int m, const int n, const int k, double *A, double *B, double *C);
extern void     mm_CTN_mmBLAS_d            (const int m, const int n, const int k, double *A, double *B, double *C);

struct S_TEST {
	// Counters
	unsigned int Ntest, Npass, Nwarnings;
};
extern struct S_TEST TestDB;

#endif // DPG__test_h__INCLUDED
