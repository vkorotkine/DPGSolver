// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__functions_h__INCLUDED
#define DPG__functions_h__INCLUDED

#include "database.h"
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
extern void     initialize_ELEMENTs         (void);
extern void     finalize_ELEMENTs           (void);
extern void     gmsh_reader                 (void);
extern void     setup_connectivity          (void);
extern void     setup_periodic              (void);
extern void       find_periodic_connections (unsigned int *Pve, unsigned int *pvePointer, const unsigned int VeMax);
extern void   setup_operators               (void);
extern void     cubature_TP                 (double **rst, double **w, unsigned int **symms, unsigned int *Nn,
                                             unsigned int *Ns, const unsigned int return_w, const unsigned int P,
                                             const unsigned int d, const char *NodeType);
extern void     cubature_TRI                (double **rst, double **w, unsigned int **symms, unsigned int *Nn,
                                             unsigned int *Ns, const unsigned int return_w, const unsigned int P,
                                             const unsigned int d, const char *NodeType);
extern void     cubature_TET                (double **rst, double **w, unsigned int **symms, unsigned int *Nn,
                                             unsigned int *Ns, const unsigned int return_w, const unsigned int P,
                                             const unsigned int d, const char *NodeType);
extern void     cubature_PYR                (double **rst, double **w, unsigned int **symms, unsigned int *Nn,
                                             unsigned int *Ns, const unsigned int return_w, const unsigned int P,
                                             const unsigned int d, const char *NodeType);
extern double     *basis_TP                 (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern double     *basis_SI                 (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern double     *basis_PYR                (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern double     **grad_basis_TP           (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern double     **grad_basis_SI           (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern double     **grad_basis_PYR          (const unsigned int P, const double *rst, const unsigned int Nn,
                                             unsigned int *NbfOut, const unsigned int d);
extern void       rst_to_abc_SI             (const unsigned int Nn, const unsigned int d, const double *rst,
                                             double *a, double *b, double *c);
extern void       rst_to_abc_PYR            (const unsigned int Nn, const unsigned int d, const double *rst,
                                             double *a, double *b, double *c);
extern double     jacobiP                   (const double x, const double alpha, const double beta, const int N);
extern double     grad_jacobiP              (const double x, const double alpha, const double beta, const int N);
extern void   setup_structures              (void);
extern void   setup_geometry                (void);
extern void     vertices_to_exact_geom      (void);
extern void     setup_ToBeCurved            (struct S_VOLUME *VOLUME);
extern void     setup_straight              (struct S_VOLUME *VOLUME);
extern void     setup_FACET_XYZ             (struct S_FACET *FACET);
extern void     setup_geom_factors          (struct S_VOLUME *VOLUME);
extern void     setup_normals               (struct S_FACET *FACET);

// Solver
extern void   compute_exact_solution (const unsigned int Nn, double *XYZ, double *UEx, double *sEx,
                                      const unsigned int solved);
extern void   initialize_test_case   (const unsigned int adapt_update_MAX);
extern void   update_VOLUME_Ops      (void);
extern void   update_VOLUME_finalize (void);
extern void   solver_explicit        (void);
extern void     explicit_VOLUME_info (void);
extern void     explicit_FACET_info  (void);
extern void       get_facet_ordering (const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
                                      const unsigned int Ns, const unsigned int Nn, const unsigned int *symms,
                                      const double *rst, unsigned int *nOrd);
extern double finalize_RHS           (void);

// Fluxes
extern void flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
                           const unsigned int Neq);
extern void flux_LF       (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);
extern void flux_ROE      (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);

// Boundary conditions
extern void boundary_Riemann  (const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut,
                               double *WB, double *nL);
extern void boundary_SlipWall (const unsigned int Nn, const unsigned int Nel, double *WL, double *WB, double *nL);

// hp adaptation
extern void         get_PS_range         (unsigned int *PSMin, unsigned int *PSMax);
extern void         get_Pb_range         (const unsigned int P, unsigned int *PbMin, unsigned int *PbMax);
extern void         get_vh_range         (const struct S_VOLUME *VOLUME, unsigned int *vhMin, unsigned int *vhMax);
extern void         get_fh_range         (const struct S_VOLUME *VOLUME, const unsigned int f, unsigned int *fhMin,
                                          unsigned int *fhMax);
extern unsigned int get_VOLUMEc_type     (const unsigned int VType, const unsigned int vh);
extern unsigned int get_IndEhref         (const unsigned int VType, const unsigned int vh);
extern void         mesh_update          (void);
extern void         mesh_to_level        (const unsigned int level);
extern void         adapt_hp             (void);
extern void           update_VOLUME_hp   (void);
extern void           update_FACET_hp    (void);
extern void           update_Vgrp        (void);
extern void           update_VOLUME_list (void);

// Sum Factorization
extern void   get_sf_parameters(const unsigned int NIn0, const unsigned int NOut0, double *OP0,
                                const unsigned int NIn1, const unsigned int NOut1, double *OP1,
                                unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                                const unsigned int d, const unsigned int dim1, const unsigned int Eclass);
extern void get_sf_parametersV(const unsigned int NIn0, const unsigned int NOut0, double **OP0, 
                               const unsigned int NIn1, const unsigned int NOut1, double **OP1, 
                               unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                               const unsigned int d, const unsigned int vh, const unsigned int Eclass);
extern void get_sf_parametersF(const unsigned int NIn0, const unsigned int NOut0, double **OP0, 
                               const unsigned int NIn1, const unsigned int NOut1, double **OP1, 
                               unsigned int NIn_SF[3], unsigned int NOut_SF[3], double *OP_SF[3],
                               const unsigned int d, const unsigned int Vf, const unsigned int Eclass);
extern void   sf_operate_d   (const unsigned int NOut, const unsigned int NCols, const unsigned int NIn,
                              const unsigned int BRowMaxIn, double *OP, double *Input, double *Output);
extern void   sf_swap_d      (double *Input, const unsigned int NRows, const unsigned int NCols,
                              const unsigned int iBound, const unsigned int jBound, const unsigned int kBound,
                              const unsigned int iStep, const unsigned int jStep, const unsigned int kStep);
extern void   sf_apply_d     (double *Input, double *Output, const unsigned int NIn[3], const unsigned int NOut[3],
                              const unsigned int NCols, double *OP[3], const unsigned int Diag[3],
                              const unsigned int d);
extern double *sf_assemble_d (const unsigned int NIn[3], const unsigned int NOut[3], const unsigned int d,
                              double *BOP[3]);

// Postprocessor

// Error computation
extern void compute_errors        (const struct S_VOLUME *VOLUME, double *L2Error2, double *Vol, unsigned int *DOF,
                                   const unsigned int solved);
extern void compute_errors_global (void);

// Plotting
extern void output_to_paraview    (const char *OutputType);
extern void plotting_element_info (double **rst, unsigned int **connect, unsigned int **types, unsigned int *Nn,
                                   unsigned int *NE, const unsigned int P, const unsigned int typeIn);

// Matrix Functions
extern double *diag_d     (const double *x, const unsigned int N);
extern double *identity_d (const unsigned int N);
extern double *inverse_d  (const unsigned int N, const unsigned int NRHS, const double *A, const double *b);
extern double *mm_Alloc_d (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                           const int m, const int n, const int k, const double alpha, const double *A, const double *B);
extern void   mm_d        (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb,
                           const int m, const int n, const int k, const double alpha, const double *A, const double *B,
                           double *C);
extern void   mm_CTN_d    (const int m, const int n, const int k, double *A, double *B, double *C);
extern void   mm_CTN_CSR_d     (const int m, const int n, const int k, const struct S_OpCSR *A, double *B, double *C);
extern void   convert_to_CSR_d (const unsigned int NRows, const unsigned int NCols, const double *Input,
                                struct S_OpCSR **Output);

// Math Functions
extern unsigned int factorial_ull (const unsigned int n);
extern unsigned int gamma_ull     (const unsigned int n);
extern double       gamma_d       (const double x);

// Struct related functions
extern struct S_ELEMENT *New_ELEMENT        (void);
extern struct S_VOLUME  *New_VOLUME         (void);
extern struct S_FACET   *New_FACET          (void);

extern unsigned int     is_ELEMENT_present  (const unsigned int type);
extern unsigned int     get_Eclass          (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_type   (const unsigned int type);
extern struct S_ELEMENT *get_ELEMENT_Eclass (const unsigned int type, const unsigned int IndEclass);
extern struct S_ELEMENT *get_ELEMENT_FACET  (const unsigned int type, const unsigned int IndEclass);
extern unsigned int     get_IndFType        (const unsigned int Eclass, const unsigned int f);

// Variable related functions
extern void convert_variables (double *VarIn, double *VarOut, const unsigned int dIn, const unsigned int dOut,
                               const unsigned int Nn, const unsigned int Nel, const char TypeIn, const char TypeOut);

// Memory Management
extern void memory_free          (void);
extern void memory_constructors  (void);
extern void memory_destructor_E  (struct S_ELEMENT *ELEMENT);
extern void memory_destructor_V  (struct S_VOLUME *VOLUME);
extern void memory_destructor_F  (struct S_FACET *FACET);
extern void memory_free_children (void);


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
	extern double       array_norm_diff_d  (const unsigned int LenA, const double *A, const double *B,
	                                        const char *NormType);

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
	extern void array_free4_d  (unsigned int iMax, unsigned int jMax, unsigned int kMax, double ****A);
	extern void array_free5_d  (unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax, double *****A);

	extern void array_free1_CSR_d (struct S_OpCSR *A);
	extern void array_free4_CSR_d (unsigned int iMax, unsigned int jMax, unsigned int kMax, struct S_OpCSR ****A);
	extern void array_free5_CSR_d (unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax,
	                               struct S_OpCSR *****A);

#endif // DPG__functions_h__INCLUDED
