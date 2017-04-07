// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_viscous_h__INCLUDED
#define DPG__fluxes_viscous_h__INCLUDED

extern void flux_viscous (const unsigned int Nn, const unsigned int Nel, const double *const W,
                          const double *const *const Q, double *const F);

#endif // DPG__fluxes_viscous_h__INCLUDED
