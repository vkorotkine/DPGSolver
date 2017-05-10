// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_h__INCLUDED
#define DPG__fluxes_inviscid_h__INCLUDED

#include "fluxes_structs.h"

extern void flux_inviscid (struct S_FLUX *const FLUXDATA);
extern void flux_LF       (const unsigned int Nn, const unsigned int Nel, const double *const WL,
                           const double *const WR, double *const nFluxNum, const double *const nL, const unsigned int d,
                           const unsigned int Neq);
extern void flux_Roe      (const unsigned int Nn, const unsigned int Nel, const double *const WL,
                           const double *const WR, double *const nFluxNum, const double *const nL, const unsigned int d,
                           const unsigned int Neq);

#endif // DPG__fluxes_inviscid_h__INCLUDED
