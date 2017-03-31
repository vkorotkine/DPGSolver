// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_c_h__INCLUDED
#define DPG__solver_functions_c_h__INCLUDED

#include <complex.h>

#include "solver_functions.h"

extern void compute_W_vI_c       (struct S_VDATA *VDATA, double complex *W_vI);
extern void convert_between_rp_c (const unsigned int Nn, const unsigned int Nrc, const double *C, double complex *Ap,
                                  double complex *Ar, const char *conv_type);

extern void finalize_VOLUME_Inviscid_Weak_c (const unsigned int Nrc, const double complex *Ar_vI, double complex *RLHS,
                                             const char *term_type, struct S_VDATA *VDATA);

#endif // DPG__solver_functions_c_h__INCLUDED
