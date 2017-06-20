// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__array_free_h__INCLUDED
#define DPG__array_free_h__INCLUDED

#include <stddef.h>
#include <complex.h>

#include "S_OpCSR.h"
#include "matrix_structs.h"

extern void free_NULL (void *A);

extern void array_free2_c     (size_t const iMax, char           **A);
extern void array_free2_ui    (size_t const iMax, unsigned int   **A);
extern void array_free2_i     (size_t const iMax, int            **A);
extern void array_free2_l     (size_t const iMax, long           **A);
extern void array_free2_ll    (size_t const iMax, long long      **A);
extern void array_free2_f     (size_t const iMax, float          **A);
extern void array_free2_d     (size_t const iMax, double         **A);
extern void array_free2_cmplx (size_t const iMax, double complex **A);
extern void array_free2_ld    (size_t const iMax, long double    **A);
extern void array_free3_c     (size_t const iMax, size_t const jMax, char         ***A);
extern void array_free3_ui    (size_t const iMax, size_t const jMax, unsigned int ***A);
extern void array_free3_i     (size_t const iMax, size_t const jMax, int          ***A);
extern void array_free3_l     (size_t const iMax, size_t const jMax, long         ***A);
extern void array_free3_ll    (size_t const iMax, size_t const jMax, long long    ***A);
extern void array_free3_f     (size_t const iMax, size_t const jMax, float        ***A);
extern void array_free3_d     (size_t const iMax, size_t const jMax, double       ***A);
extern void array_free3_ld    (size_t const iMax, size_t const jMax, long double  ***A);
extern void array_free4_ui    (size_t const iMax, size_t const jMax, size_t const kMax, unsigned int ****A);
extern void array_free4_d     (size_t const iMax, size_t const jMax, size_t const kMax, double       ****A);
extern void array_free5_d     (size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, double *****A);

extern void array_free1_CSR_d (struct S_OpCSR *A);
extern void array_free4_CSR_d (size_t const iMax, size_t const jMax, size_t const kMax, struct S_OpCSR ****A);
extern void array_free5_CSR_d (size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, struct S_OpCSR *****A);

extern void matrix_free  (struct S_MATRIX *A);
extern void matrix_free2 (size_t const iMax, struct S_MATRIX **A);
extern void matrix_free3 (size_t const N0, size_t const N1, struct S_MATRIX ***A);
extern void matrix_free4 (size_t const N0, size_t const N1, size_t const N2, struct S_MATRIX ****A);
extern void matrix_free5 (size_t const N0, size_t const N1, size_t const N2, size_t const N3, struct S_MATRIX *****A);

extern void multiarray_free (struct S_MULTI_ARRAY *A);

#endif // DPG__array_free_h__INCLUDED
