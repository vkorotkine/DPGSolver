// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_structs_h__INCLUDED
#define DPG__fluxes_structs_h__INCLUDED

#include <complex.h>

struct S_FLUX {
	unsigned int d, Nn, Nel;

	double const *W,
	             *const *Q;
	double *F, *Fr, *dFdW, *dFrdW, **dFdQ, **dFrdQ;

	// Only used for verification.
	double complex const *W_c,
	                     *const *Q_c;
	double complex *F_c, *Fr_c;
};

#endif // DPG__fluxes_structs_h__INCLUDED
