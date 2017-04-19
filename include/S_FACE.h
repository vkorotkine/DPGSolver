// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_FACE_h__INCLUDED
#define DPG__S_FACE_h__INCLUDED

#include <complex.h>

struct S_FACE {
	// Structures
	unsigned int P, type, VfIn, VfOut, indexg, BC, IndOrdInOut, IndOrdOutIn, level, update, adapt_type;

	// Geometry
	unsigned int curved, typeInt;
	double       *XYZ_fI, *XYZ_fS, *n_fI, *n_fS, *detJF_fI, *detJF_fS, *detJVIn_fI, *detJVOut_fI;

	// Solving
	unsigned int Boundary;
	double       *RHSIn, *RHSOut, *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut,
	             **QhatL, **QhatR, **Qhat_WhatLL, **Qhat_WhatRL, **Qhat_WhatLR, **Qhat_WhatRR;

	// Poisson
	double **qhatIn, **qhatOut, **qhat_uhatInIn, **qhat_uhatOutIn, **qhat_uhatInOut, **qhat_uhatOutOut;

	// Linearization testing
	double complex *RHSIn_c, *RHSOut_c, **qhatIn_c, **qhatOut_c, **QhatL_c, **QhatR_c;

	// structs
	struct S_VOLUME *VIn, *VOut;
	struct S_FACE  *next, *grpnext, *child0, *parent;
};

#endif // DPG__S_FACE_h__INCLUDED
