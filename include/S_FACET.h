// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__S_FACET_h__INCLUDED
#define DPG__S_FACET_h__INCLUDED

#include <complex.h>

struct S_FACET {
	// Structures
	unsigned int P, type, VfIn, VfOut, indexg, BC, IndOrdInOut, IndOrdOutIn, level, update, adapt_type;

	// Geometry
	unsigned int curved, typeInt;
	double       *XYZ_fI, *XYZ_fS, *n_fI, *n_fS, *detJF_fI, *detJF_fS;

	// Solving
	unsigned int Boundary;
	double       *RHSIn, *RHSOut, *LHSInIn, *LHSOutIn, *LHSInOut, *LHSOutOut;

	// Linearization testing
	double complex *RHSIn_c, *RHSOut_c;

	// structs
	struct S_VOLUME *VIn, *VOut;
	struct S_FACET  *next, *grpnext, *child0, *parent;
};

#endif // DPG__S_FACET_h__INCLUDED
