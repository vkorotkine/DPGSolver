// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__jacobian_fluxes_inviscid_h__INCLUDED
#define DPG__jacobian_fluxes_inviscid_h__INCLUDED

extern void jacobian_flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *dFdW,
                                    const unsigned int d, const unsigned int Neq);

#endif // DPG__jacobian_fluxes_inviscid_h__INCLUDED
