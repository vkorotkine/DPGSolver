// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_viscous_h__INCLUDED
#define DPG__fluxes_viscous_h__INCLUDED

extern void flux_viscous (unsigned int const Nn, unsigned int const Nel, double const *const W,
                          double const *const *const Q, double *const F);

#endif // DPG__fluxes_viscous_h__INCLUDED
