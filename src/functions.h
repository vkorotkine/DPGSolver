#ifndef DPG__functions_h__INCLUDED
#define DPG__functions_h__INCLUDED

#include "mkl.h"

/*
 *	Purpose:
 *		Set function prototypes.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

// Preprocessor
extern void   initialization                (int nargc, char **argv);
extern void   setup_parameters              (void);
extern void   setup_mesh                    (void);
extern void     gmsh_reader                 (void);
extern void     setup_connectivity          (void);
extern void     setup_periodic              (void);
extern void       find_periodic_connections (unsigned int *Pve, unsigned int *pvePointer, const unsigned int VeMax);
extern void   setup_operators               (void);
extern void     cubature_TP                 (double **xir, double **W, unsigned int **Con, unsigned int *Nn,
                                             const unsigned int *ToReturn, const unsigned int P, const unsigned int d,
											 const char *NodeType);
extern double     *basis_TP                 (const int P, const double *xir, const int Nn, const int d);
extern double     **grad_basis_TP           (const int P, const double *xir, const int Nn, const int d);
extern double     jacobiP                   (const double x, const double alpha, const double beta, const int N);
extern double     grad_jacobiP              (const double x, const double alpha, const double beta, const int N);
extern void   setup_structures              (void);
extern void   setup_geometry                (void);
extern void     vertices_to_exact_geom      (void);

// Matrix Functions
extern double *identity_d (const int N);
extern double *inverse_d  (int N, int NRHS, double *A, double *b);
extern double *mm_Alloc_d (const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m, const int n, const int k,
                           const double alpha, const double *A, const double *B);
extern void   mm_d        (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const int m,
                           const int n, const int k, const double alpha, const double *A, const double *B, double *C);
extern void   mm_CTN_d    (const int m, const int n, const int k, const double *A, const double *B, double *C);

//extern void   mm_CTN_d    (const int m, const int n, const int k, const double *A, const double *B, double *C,
//                           const int useBLAS);

// Math Functions
extern int    factorial_i (const int n);
extern int    gamma_i     (const int n);
extern double gamma_d     (const double x);

// Struct related functions
extern struct S_ELEMENT *New_ELEMENT        (void);
extern struct S_VOLUME  *New_VOLUME         (void);

extern int              is_ELEMENT_present  (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_type   (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_Eclass (const unsigned int Eclass, const unsigned int Esubclass);


// Memory Management
extern void memory_free         (void);
extern void memory_constructors (void);
extern void memory_destructor_E (struct S_ELEMENT *ELEMENT);


// Array Processing
	// Sorting
	extern void array_sort_ui        (unsigned int NRows, unsigned int NCols, unsigned int *A, unsigned int *Indices,
	                                  const char ordering, const char trans);
	extern void array_sort_i         (unsigned int NRows, unsigned int NCols, int *A, unsigned int *Indices,
	                                  const char ordering, const char trans);
	extern void array_sort_d         (unsigned int NRows, unsigned int NCols, double *A, unsigned int *Indices,
	                                  const char ordering, const char trans);
	extern void array_find_indexo_ui (const unsigned int LenA, const unsigned int *A, const unsigned int val, unsigned int *IdxF, unsigned int *LenF);

	// Norms
	extern double array_norm_d (int LenA, double *A, char *NormType);

	// Swapping
	extern void array_swap_ui (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
	                           const unsigned int stepIn);
	extern void array_swap_i  (register int *arr1, register int *arr2, const int NIn, const int stepIn);
	extern void array_swap_d  (register double *arr1, register double *arr2, const int NIn, const int stepIn);

	// Printing
//	extern void array_print    (const unsigned int m, const unsigned int n, void *A, char *type);
	extern void array_print_ui (const unsigned int m, const unsigned int n, const unsigned int *A, const char layout);
	extern void array_print_i  (const unsigned int m, const unsigned int n, int *A, const char layout);
	extern void array_print_l  (const unsigned int m, const unsigned int n, long *A, const char layout);
	extern void array_print_ll (const unsigned int m, const unsigned int n, long long *A, const char layout);
	extern void array_print_f  (const unsigned int m, const unsigned int n, float *A, const char layout);
	extern void array_print_d  (const unsigned int m, const unsigned int n, double *A, const char layout);
	extern void array_print_ld (const unsigned int m, const unsigned int n, long double *A, const char layout);

	// Memory Management
	extern void array_free2_c  (int iMax, char **A);
	extern void array_free2_i  (int iMax, int **A);
	extern void array_free2_l  (int iMax, long **A);
	extern void array_free2_ll (int iMax, long long **A);
	extern void array_free2_f  (int iMax, float **A);
	extern void array_free2_d  (int iMax, double **A);
	extern void array_free2_ld (int iMax, long double **A);
	extern void array_free3_c  (int iMax, int jMax, char ***A);
	extern void array_free3_i  (int iMax, int jMax, int ***A);
	extern void array_free3_l  (int iMax, int jMax, long ***A);
	extern void array_free3_ll (int iMax, int jMax, long long ***A);
	extern void array_free3_f  (int iMax, int jMax, float ***A);
	extern void array_free3_d  (int iMax, int jMax, double ***A);
	extern void array_free3_ld (int iMax, int jMax, long double ***A);

// Testing
	extern void test_print                (const unsigned int pass);
	extern void test_speed_mm_d           (void);
	extern void test_imp_array_find_index (void);
	extern void test_imp_array_norm       (void);
	extern void test_imp_array_sort       (void);

#endif // DPG_functions_h__INCLUDED
