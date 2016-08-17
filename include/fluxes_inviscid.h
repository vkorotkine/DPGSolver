// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__fluxes_inviscid_h__INCLUDED
#define DPG__fluxes_inviscid_h__INCLUDED

extern void flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
                           const unsigned int Neq);
extern void flux_LF       (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);
extern void flux_Roe      (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);

#endif // DPG__fluxes_inviscid_h__INCLUDED
