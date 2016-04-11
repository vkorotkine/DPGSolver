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

extern void test_print                         (const unsigned int pass);
extern void test_speed_mm_d                    (void);
extern void test_imp_array_find_index          (void);
extern void test_imp_array_norm                (void);
extern void test_imp_array_sort                (void);
extern void test_imp_array_swap                (void);
extern void test_imp_math_factorial            (void);
extern void test_imp_math_gamma                (void);
extern void test_imp_matrix_identity           (void);
extern void test_imp_matrix_inverse            (void);
extern void test_imp_matrix_mm                 (void);
extern void test_imp_find_periodic_connections (void);
extern void test_imp_cubature_TP               (void);
extern void test_imp_basis_TP                  (void);

struct S_TEST {
	// Counters
	unsigned int Ntest, Npass, Nwarnings;
};
extern struct S_TEST TestDB;

#endif // DPG__test_h__INCLUDED
