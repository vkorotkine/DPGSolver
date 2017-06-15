// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__jacobian_fluxes_inviscid_h__INCLUDED
#define DPG__jacobian_fluxes_inviscid_h__INCLUDED

#include "fluxes_structs.h"

extern void jacobian_flux_inviscid     (struct S_FLUX *const FLUXDATA);
extern void jacobian_flux_num_inviscid (struct S_NUMERICALFLUX *const NUMFLUXDATA);

extern void jacobian_flux_inviscid_M (struct S_FLUX_M *const FLUXDATA_M);

#endif // DPG__jacobian_fluxes_inviscid_h__INCLUDED
