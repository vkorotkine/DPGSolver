// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_c_h__INCLUDED
#define DPG__fluxes_inviscid_c_h__INCLUDED

#include <complex.h>


extern void flux_inviscid_c (const unsigned int Nn, const unsigned int Nel, double complex *W, double complex *F,
                             const unsigned int d, const unsigned int Neq);
extern void flux_LF_c       (const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WR,
                             double complex *nFluxNum, double *nL, const unsigned int d, const unsigned int Neq);
extern void flux_Roe_c      (const unsigned int Nn, const unsigned int Nel, double complex *WL, double complex *WR,
                             double complex *nFluxNum, double *nL, const unsigned int d, const unsigned int Neq);

#endif // DPG__fluxes_inviscid_c_h__INCLUDED
