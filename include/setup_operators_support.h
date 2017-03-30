// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__setup_operators_support_h__INCLUDED
#define DPG__setup_operators_support_h__INCLUDED

#include "S_ELEMENT.h"

struct S_BCOORDS {
	unsigned int Nve,
	             *NfnG2, *NfnGc, *NfnS, *NfnIs, *NfnIc, *NenG2, *NenGc;
	double       **w_fIs, **w_fIc, **BCoords_G2, **BCoords_Gc, **BCoords_S, **BCoords_Is, **BCoords_Ic;
};

extern struct S_BCOORDS *get_BCoords_dEm1 (const struct S_ELEMENT *ELEMENT, const unsigned int IndFType);
extern struct S_BCOORDS *get_BCoords_dEm2 (const struct S_ELEMENT *ELEMENT);

extern void   setup_ELEMENT_plotting      (const unsigned int EType);
extern void   setup_ELEMENT_normals       (const unsigned int EType);
extern double *get_rst_vV                 (const struct S_ELEMENT *ELEMENT);
extern void   setup_ELEMENT_VeV           (const unsigned int EType);
extern void   setup_ELEMENT_VeF           (const unsigned int EType);
extern void   setup_ELEMENT_VeE           (const unsigned int EType);
extern double get_L2_scaling              (const unsigned int EType, const unsigned int vref);
extern void   get_face_ordering           (const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
                                           const unsigned int Ns, const unsigned int Nn, const unsigned int *symms,
                                           const double *rst, unsigned int *nOrd);
extern void   setup_ELEMENT_FACE_ordering (const unsigned int FType);
extern void   compute_ELEMENT_Volume      (const unsigned int EType);

#endif // DPG__setup_operators_support_h__INCLUDED
