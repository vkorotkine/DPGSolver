// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__array_free_h__INCLUDED
#define DPG__array_free_h__INCLUDED

#include <complex.h>

#include "S_OpCSR.h"


extern void array_free2_c     (unsigned int iMax, char           **A);
extern void array_free2_ui    (unsigned int iMax, unsigned int   **A);
extern void array_free2_i     (unsigned int iMax, int            **A);
extern void array_free2_l     (unsigned int iMax, long           **A);
extern void array_free2_ll    (unsigned int iMax, long long      **A);
extern void array_free2_f     (unsigned int iMax, float          **A);
extern void array_free2_d     (unsigned int iMax, double         **A);
extern void array_free2_cmplx (unsigned int iMax, double complex **A);
extern void array_free2_ld    (unsigned int iMax, long double    **A);
extern void array_free3_c     (unsigned int iMax, unsigned int jMax, char         ***A);
extern void array_free3_ui    (unsigned int iMax, unsigned int jMax, unsigned int ***A);
extern void array_free3_i     (unsigned int iMax, unsigned int jMax, int          ***A);
extern void array_free3_l     (unsigned int iMax, unsigned int jMax, long         ***A);
extern void array_free3_ll    (unsigned int iMax, unsigned int jMax, long long    ***A);
extern void array_free3_f     (unsigned int iMax, unsigned int jMax, float        ***A);
extern void array_free3_d     (unsigned int iMax, unsigned int jMax, double       ***A);
extern void array_free3_ld    (unsigned int iMax, unsigned int jMax, long double  ***A);
extern void array_free4_d     (unsigned int iMax, unsigned int jMax, unsigned int kMax, double ****A);
extern void array_free5_d     (unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax, double *****A);

extern void array_free1_CSR_d (struct S_OpCSR *A);
extern void array_free4_CSR_d (unsigned int iMax, unsigned int jMax, unsigned int kMax, struct S_OpCSR ****A);
extern void array_free5_CSR_d (unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax, struct S_OpCSR *****A);

#endif // DPG__array_free_h__INCLUDED
