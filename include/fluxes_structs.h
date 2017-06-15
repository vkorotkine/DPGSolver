// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_structs_h__INCLUDED
#define DPG__fluxes_structs_h__INCLUDED

#include <complex.h>
#include "matrix_structs.h"

struct S_FLUX {
	unsigned int d, Nn, Nel, PDE_index;

	double const *W, *XYZ,
	             *const *Q;
	double *F, *Fr, *dFdW, *dFrdW, **dFdQ, **dFrdQ;

	// Only used for verification.
	double complex const *W_c,
	                     *const *Q_c;
	double complex *F_c, *Fr_c;
};

struct S_NUMERICALFLUX {
	unsigned int d, Nn, Nel, NumFluxInviscid_index, NumFluxViscous_index;

	double const *nL, *XYZ,
	             *WL, *WR;
	double       *nFluxNum, *dnFluxNumdWL, *dnFluxNumdWR, **dnFluxNumdQL, **dnFluxNumdQR,
	             **nSolNum, **dnSolNumdWL, **dnSolNumdWR;

	// Only used for verification.
	double complex const *WL_c, *WR_c;
	double complex       *nFluxNum_c, **nSolNum_c;
};

struct S_FLUX_M {
	unsigned int d;

	struct S_MATRIX const *W, *XYZ;
	struct S_MATRIX *F, *Fr, *dFdW, *dFrdW;
};

#endif // DPG__fluxes_structs_h__INCLUDED
