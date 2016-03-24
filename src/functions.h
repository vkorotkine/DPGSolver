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

#ifndef DPG__functions_h__INCLUDED
#define DPG__functions_h__INCLUDED

// Preprocessor
extern void   initialization                (int nargc, char **argv);
extern void   setup_parameters              (void);
extern void   setup_mesh                    (void);
extern void     gmsh_reader                 (void);
extern void     setup_connectivity          (void);
extern void     setup_periodic              (void);
extern void       find_periodic_connections (int *Pve, int *pvePointer, int VeMax);
extern void   setup_operators               (void);
extern void     cubature_TP                 (double **xir, double **W, int **Con, int *Nn, int *ToReturn, int P, int d, char *NodeType);
extern double     *basis_TP                 (const int P, const double *xir, const int Nn, const int d);
extern double     **grad_basis_TP           (const int P, const double *xir, const int Nn, const int d);
extern double     jacobiP                   (const double x, const double alpha, const double beta, const int N);
extern double     grad_jacobiP              (const double x, const double alpha, const double beta, const int N);

// Math Functions
extern int    factorial_i (const int n);
extern int    gamma_i     (const int n);
extern double gamma_d     (const double x);

// Extern structs
extern struct S_ELEMENT *New_ELEMENT (void);


// Memory Management
extern void memory_free         (void);
extern void memory_constructors (void);
extern void memory_destructor_E (struct S_ELEMENT *ELEMENT);


// Array Processing
	// Sorting
	extern void array_sort_i        (int NRows, int NCols, int    *A, int *Indices, char ordering, char trans);
	extern void array_sort_d        (int NRows, int NCols, double *A, int *Indices, char ordering, char trans);
	extern void array_find_indexo_i (int LenA, int *A, int val, int *IdxF, int *LenF);

	// Norms
	extern double array_norm_d (int LenA, double *A, char *NormType);

	// Printing
//	extern void array_print    (int m, int n, void *A, char *type);
	extern void array_print_i  (int m, int n, int *A);
	extern void array_print_l  (int m, int n, long *A);
	extern void array_print_ll (int m, int n, long long *A);
	extern void array_print_f  (int m, int n, float *A);
	extern void array_print_d  (int m, int n, double *A);
	extern void array_print_ld (int m, int n, long double *A);

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


#endif // DPG_functions_h__INCLUDED
