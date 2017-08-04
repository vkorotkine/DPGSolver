// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_VOLUME_h__INCLUDED
#define DPG__S_VOLUME_h__INCLUDED
/// \file

#include <complex.h>

#include "matrix_structs.h"
#include "S_ELEMENT.h"

/**\{ \name Alternate notation
 *  Provide an alternate notation to be deleted after eventual code refactoring.
 *  \todo Fix the notation.
 */
#define Volume S_VOLUME
///\}

/** \brief Struct holding data related to the base Volume. */
struct S_VOLUME {
	// Structures
	unsigned int indexl, indexg, P, type, Eclass, update, curved, level,
	             *NsubF, *neigh, *neigh_f;
	double *XYZ_vV, *XYZ_vVc;

	struct S_ELEMENT*const element;

	// Geometry
	unsigned int NvnG, *VeInd, *VeInfo, **BC;
	double *XYZ_S, *XYZ, *detJV_vI, *C_vC, *C_vI, **C_vf;

	struct S_MULTI_ARRAY *XYZ_MA, *detJV_vI_MA, *C_vI_MA;

	// Initialization
	unsigned int NvnS;
	double *What, **QhatV, **Qhat, *RES, **QhatV_What, **Qhat_What;

	struct S_MULTI_ARRAY *What_MA;

	// Solving
	unsigned int   IndA, nnz_d, nnz_o;
	double         *RHS, *LHS, **LHSQ, *wdetJV_vI, *MInv, *MInv_diag;

	struct S_MULTI_ARRAY *RHS_MA, *LHS_L_MA, *LHS_R_MA;

	// Linearization testing
	double complex *RHS_c, *What_c, **Qhat_c, **QhatV_c, **qhat_c;

	// hp adaptivity
	unsigned int refine_current, Vadapt, adapt_type, PNew, hrefine_type;
	double       *XYZ_vVP2;
//	double       minRES, maxRES;

	// structs
	struct S_VOLUME *next, *grpnext, *child0, *parent;
	struct S_FACE  **FACE;
};

void set_element (struct S_VOLUME*const volume, struct S_ELEMENT*const element_head);

#endif // DPG__S_VOLUME_h__INCLUDED
