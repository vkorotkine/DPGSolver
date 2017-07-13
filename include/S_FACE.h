// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_FACE_h__INCLUDED
#define DPG__S_FACE_h__INCLUDED

#include <complex.h>

struct S_FACE {
	// Structures
	unsigned int P, type, VfL, VfR, indexg, BC, IndOrdLR, IndOrdRL, level, update, adapt_type;

	// Geometry
	unsigned int curved, typeInt;
	double       *XYZ_fI, *XYZ_fS, *n_fI, *n_fS, *detJF_fI, *detJF_fS, *detJVIn_fI, *detJVOut_fI;

	// Solving
	char         CDG2_side;
	unsigned int Boundary;
	double       *RHSL, *RHSR, *LHSLL, *LHSRL, *LHSLR, *LHSRR,
	             **QhatL, **QhatR, **QhatL_WhatL, **QhatL_WhatR, **QhatR_WhatL, **QhatR_WhatR;

	// Linearization testing
	double complex *RHSL_c, *RHSR_c, **QhatL_c, **QhatR_c;

	// structs
	struct S_VOLUME *VL, *VR;
	struct S_FACE  *next, *grpnext, *child0, *parent;
};

#endif // DPG__S_FACE_h__INCLUDED
