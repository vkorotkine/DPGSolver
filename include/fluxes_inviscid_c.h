// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_c_h__INCLUDED
#define DPG__fluxes_inviscid_c_h__INCLUDED

#include <complex.h>

#include "fluxes_structs.h"

extern void flux_inviscid_c (struct S_FLUX *const FLUX);
extern void flux_LF_c       (const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
                             const double complex *const WR, double complex *const nFluxNum, const double *const nL,
                             const unsigned int d, const unsigned int Neq);
extern void flux_Roe_c      (const unsigned int Nn, const unsigned int Nel, const double complex *const WL,
                             const double complex *const WR, double complex *const nFluxNum, const double *const nL,
                             const unsigned int d, const unsigned int Neq);

#endif // DPG__fluxes_inviscid_c_h__INCLUDED
