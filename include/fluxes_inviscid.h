// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_h__INCLUDED
#define DPG__fluxes_inviscid_h__INCLUDED

#include "fluxes_structs.h"

extern void flux_inviscid     (struct S_FLUX *const FLUXDATA);
extern void flux_num_inviscid (struct S_NUMERICALFLUX *const NUMFLUXDATA);
extern void flux_Euler        (struct S_FLUX *const FLUXDATA);

extern void flux_inviscid_M (struct S_FLUX_M *const FLUXDATA_M);

#endif // DPG__fluxes_inviscid_h__INCLUDED
