#ifndef DPG__functions_h__INCLUDED
#define DPG__functions_h__INCLUDED

#include "mkl.h"

/*
 *	Purpose:
 *		Set function prototypes.
 *
 *	Comments:
 *		ToBeDeleted: Ensure that lines are not longer than 120 characters.
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
extern void     cubature_TP                 (double **rst, double **w_vec, unsigned int **Con, unsigned int *Nn,
                                             const unsigned int *ToReturn, const unsigned int P, const unsigned int d,
											 const char *NodeType);
extern double     *basis_TP                 (const unsigned int P, const double *rst, const unsigned int Nn,
                                             const unsigned int d);
extern double     **grad_basis_TP           (const unsigned int P, const double *rst, const unsigned int Nn,
                                             const unsigned int d);
extern double     jacobiP                   (const double x, const double alpha, const double beta, const int N);
extern double     grad_jacobiP              (const double x, const double alpha, const double beta, const int N);
extern void   setup_structures              (void);
extern void   setup_geometry                (void);
extern void     vertices_to_exact_geom      (void);

// Matrix Functions
extern double *diag_d     (const double *x, const unsigned int N);
extern double *identity_d (const unsigned int N);
extern double *inverse_d  (const unsigned int N, const unsigned int NRHS, const double *A, const double *b);
extern double *mm_Alloc_d (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                           const int m, const int n, const int k, const double alpha, const double *A, const double *B);
extern void   mm_d        (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                           const int m, const int n, const int k, const double alpha, const double *A, const double *B,
                           double *C);
extern void   mm_CTN_d    (const int m, const int n, const int k, const double *A, const double *B, double *C);

//extern void   mm_CTN_d    (const int m, const int n, const int k, const double *A, const double *B, double *C,
//                           const int useBLAS);

// Math Functions
extern unsigned int factorial_ui (const unsigned int n);
extern unsigned int gamma_ui     (const unsigned int n);
extern double       gamma_d      (const double x);

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
	extern void array_find_indexo_ui (const unsigned int LenA, const unsigned int *A, const unsigned int val,
	                                  unsigned int *IdxF, unsigned int *LenF);

	// Norms
	extern unsigned int array_norm_ui      (const unsigned int LenA, const unsigned int *A, const char *NormType);
	extern double       array_norm_d       (const unsigned int LenA, const double *A, const char *NormType);
	extern unsigned int array_norm_diff_ui (const unsigned int LenA, const unsigned int *A, const unsigned int *B,
	                                        const char *NormType);
	extern double       array_norm_diff_d  (const unsigned int LenA, const double *A, const double *B, const char *NormType);

	// Swapping
	extern void array_swap_ui (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
	                           const unsigned int stepIn);
//	extern void array_swap_i  (register int *arr1, register int *arr2, const int NIn, const int stepIn);
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
	extern void array_free2_c  (unsigned int iMax, char **A);
	extern void array_free2_ui (unsigned int iMax, unsigned int **A);
	extern void array_free2_i  (unsigned int iMax, int **A);
	extern void array_free2_l  (unsigned int iMax, long **A);
	extern void array_free2_ll (unsigned int iMax, long long **A);
	extern void array_free2_f  (unsigned int iMax, float **A);
	extern void array_free2_d  (unsigned int iMax, double **A);
	extern void array_free2_ld (unsigned int iMax, long double **A);
	extern void array_free3_c  (unsigned int iMax, unsigned int jMax, char ***A);
	extern void array_free3_ui (unsigned int iMax, unsigned int jMax, unsigned int ***A);
	extern void array_free3_i  (unsigned int iMax, unsigned int jMax, int ***A);
	extern void array_free3_l  (unsigned int iMax, unsigned int jMax, long ***A);
	extern void array_free3_ll (unsigned int iMax, unsigned int jMax, long long ***A);
	extern void array_free3_f  (unsigned int iMax, unsigned int jMax, float ***A);
	extern void array_free3_d  (unsigned int iMax, unsigned int jMax, double ***A);
	extern void array_free3_ld (unsigned int iMax, unsigned int jMax, long double ***A);

#endif // DPG__functions_h__INCLUDED
