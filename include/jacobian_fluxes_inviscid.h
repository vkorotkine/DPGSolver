// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__jacobian_fluxes_inviscid_h__INCLUDED
#define DPG__jacobian_fluxes_inviscid_h__INCLUDED

extern void jacobian_flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *dFdW,
                                    const unsigned int d, const unsigned int Neq);
extern void jacobian_flux_LF       (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *dnFdW,
                                    double *nL, const unsigned int d, const unsigned int Neq, const char side);
extern void jacobian_flux_Roe      (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *dnFdW,
                                    double *nL, const unsigned int d, const unsigned int Neq, const char side);

#endif // DPG__jacobian_fluxes_inviscid_h__INCLUDED
