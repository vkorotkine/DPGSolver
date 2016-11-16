// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_operators_support_h__INCLUDED
#define DPG__setup_operators_support_h__INCLUDED

#include "S_ELEMENT.h"

struct S_BCOORDS {
	unsigned int Nve,
	             *NfnGc, *NfnS, *NfnIs, *NfnIc, *NenGc;
	double       **w_fIs, **w_fIc, **BCoords_Gc, **BCoords_S, **BCoords_Is, **BCoords_Ic;
};


extern void             setup_ELEMENT_plotting       (const unsigned int EType);
extern void             setup_ELEMENT_normals        (const unsigned int EType);
extern double           *get_rst_vV                  (const struct S_ELEMENT *ELEMENT);
extern struct S_BCOORDS *get_BCoords_dEm1            (const struct S_ELEMENT *ELEMENT, const unsigned int IndFType);
extern struct S_BCOORDS *get_BCoords_dEm2            (const struct S_ELEMENT *ELEMENT);
extern void             setup_ELEMENT_VeV            (const unsigned int EType);
extern void             setup_ELEMENT_VeF            (const unsigned int EType);
extern void             setup_ELEMENT_VeE            (const unsigned int EType);
extern double           get_L2_scaling               (const unsigned int EType, const unsigned int vref);
extern void             setup_ELEMENT_FACE_ordering (const unsigned int FType);

#endif // DPG__setup_operators_support_h__INCLUDED
