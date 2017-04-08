// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_functions_c_h__INCLUDED
#define DPG__solver_functions_c_h__INCLUDED

#include <complex.h>

#include "solver_functions.h"

extern void coef_to_values_vI_c  (struct S_VDATA const *const VDATA, char const coef_type);
extern void convert_between_rp_c (unsigned int const Nn, unsigned int const Nrc, double const *const C,
                                  double complex *const Ap, double complex *const Ar, char const *const conv_type);

extern void finalize_VOLUME_Inviscid_Weak_c (unsigned int const Nrc, double complex const *const Ar_vI,
                                             double complex *const RLHS, char const imex_type,
                                             struct S_VDATA const *const VDATA);


extern void coef_to_values_fI_c           (struct S_FDATA const *const FDATA, char const coef_type);
extern void compute_WR_fIL_c              (struct S_FDATA const *const FDATA, double complex const *const WL_fIL,
                                           double complex *const WR_fIL);
extern void compute_numerical_flux_c      (struct S_FDATA const *const FDATA, char const imex_type);
extern void add_Jacobian_scaling_FACE_c   (struct S_FDATA const *const FDATA, char const imex_type);
extern void finalize_FACE_Inviscid_Weak_c (struct S_FDATA const *const FDATAL, struct S_FDATA const *const FDATAR,
                                           char const side, char const imex_type);

#endif // DPG__solver_functions_c_h__INCLUDED
