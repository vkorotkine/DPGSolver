// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__S_VOLUME_h__INCLUDED
#define DPG__S_VOLUME_h__INCLUDED

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, update, curved, level,
	             *NsubF, *neigh, *neigh_f;
	double *XYZ_vC;

	// Geometry
	unsigned int NvnG;
	double *XYZ_S, *XYZ, *detJV_vI, *C_vC, *C_vI, **C_vf;

	// Initialization
	unsigned int NvnS;
	double *What, *RES;

	// Solving
	double *RHS, *LHS, *wdetJV_vI, *MInv;

	// hp adaptivity
	unsigned int refine_current, Vadapt, adapt_type, PNew, hrefine_type;
//	double       minRES, maxRES;

	// structs
	struct S_VOLUME *next, *grpnext, *child0, *parent;
	struct S_FACET  **FACET;
};

#endif // DPG__S_VOLUME_h__INCLUDED
