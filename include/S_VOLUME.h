// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__S_VOLUME_h__INCLUDED
#define DPG__S_VOLUME_h__INCLUDED

#include <complex.h>

struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, update, curved, level,
	             *NsubF, *neigh, *neigh_f;
	double *XYZ_vV;

	// Geometry
	unsigned int NvnG, *VeInd, *VeInfo, **BC;
	double *XYZ_S, *XYZ, *detJV_vI, *C_vC, *C_vI, **C_vf;

	// Initialization
	unsigned int NvnS;
	double *What, *RES;

	// Solving
	unsigned int   IndA, nnz_d, nnz_o;
	double         *RHS, *LHS, *wdetJV_vI, *MInv;

	// Linearization testing
	double complex *RHS_c, *What_c, *uhat_c, **qhat_c;

	// hp adaptivity
	unsigned int refine_current, Vadapt, adapt_type, PNew, hrefine_type;
//	double       minRES, maxRES;

	// Poisson
	double *uhat, **qhat, **qhat_uhat;
	double **DxyzChiS;

	// structs
	struct S_VOLUME *next, *grpnext, *child0, *parent;
	struct S_FACE  **FACE;
};

#endif // DPG__S_VOLUME_h__INCLUDED
