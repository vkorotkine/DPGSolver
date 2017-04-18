// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_viscous_h__INCLUDED
#define DPG__fluxes_viscous_h__INCLUDED

#include <complex.h>

extern void flux_viscous_c (unsigned int const Nn, unsigned int const Nel, double complex const *const W,
                            double complex const *const *const Q, double complex *const F);

#endif // DPG__fluxes_viscous_h__INCLUDED
