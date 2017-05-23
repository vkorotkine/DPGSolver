// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__jacobian_fluxes_inviscid_h__INCLUDED
#define DPG__jacobian_fluxes_inviscid_h__INCLUDED

#include "fluxes_structs.h"

extern void jacobian_flux_inviscid     (struct S_FLUX *const FLUXDATA);
extern void jacobian_flux_num_inviscid (struct S_NUMERICALFLUX *const NUMFLUXDATA);

#endif // DPG__jacobian_fluxes_inviscid_h__INCLUDED
