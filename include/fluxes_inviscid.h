// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__fluxes_inviscid_h__INCLUDED
#define DPG__fluxes_inviscid_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"

#include "variable_functions.h"


extern void flux_inviscid (const unsigned int Nn, const unsigned int Nel, double *W, double *F, const unsigned int d,
                           const unsigned int Neq);
extern void flux_LF       (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);
extern void flux_ROE      (const unsigned int Nn, const unsigned int Nel, double *WL, double *WR, double *nFluxNum,
                           double *nL, const unsigned int d, const unsigned int Neq);

#endif // DPG__fluxes_inviscid_h__INCLUDED
