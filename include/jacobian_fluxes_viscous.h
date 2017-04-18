// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__jacobian_fluxes_viscous_h__INCLUDED
#define DPG__jacobian_fluxes_viscous_h__INCLUDED

extern void jacobian_flux_viscous (unsigned int const Nn, unsigned int const Nel, double const *const W,
	                               double const *const *const Q, double *const dFdW, double *const *const dFdQ);

#endif // DPG__jacobian_fluxes_viscous_h__INCLUDED
