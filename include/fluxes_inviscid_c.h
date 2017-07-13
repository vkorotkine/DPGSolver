// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_c_h__INCLUDED
#define DPG__fluxes_inviscid_c_h__INCLUDED

#include <complex.h>

#include "fluxes_structs.h"

extern void flux_inviscid_c     (struct S_FLUX *const FLUX);
extern void flux_num_inviscid_c (struct S_NUMERICALFLUX *const NUMFLUXDATA);

#endif // DPG__fluxes_inviscid_c_h__INCLUDED
