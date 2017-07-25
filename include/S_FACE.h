// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_FACE_h__INCLUDED
#define DPG__S_FACE_h__INCLUDED

/*
 *	Purpose:
 *		Provide a (S_)truct to hold information relating to FACEs.
 *
 *	Comments:
 *		A FACE is a d-1 dimensional ELEMENT which exists on each face of the VOLUMEs in the mesh.
 *
 *	Notation: (ToBeModified: Likely needs to be updated)
 *		Please consult Parameters.h for the definition of the constants used below.
 *		Integer division is purposely used in several definitions below.
 *		Please consult the conversion functions in S_FACE.c for equivalence between old and new notation. (ToBeDeleted)
 *
 *		P    : (P)olynomial order
 *		type : ELEMENT (type)
 *
 *		data(L/R) : (data) specific to one of the (L)eft or (R)ight VOLUMEs.
 *
 *		Ind() : (Ind)ex
 *			(fh)  : (f)ace index of potentially (h)-refined FACE as seen from the neighbouring VOLUME.
 *			        Range: [0:NFMAX*NFREFMAX)
 *			        In the case of a conforming mesh (h-refinement disabled), fh is restricted to a multiple of
 *			        NFREFMAX.
 *			(mf)  : (m)acro (f)ace number of the neighbouring VOLUME
 *			(lfh) : (l)ocal (fh) species the local h-refinement index
 *			(sfh) : (s)ub (fh) index specifies the
 *			        Range: [0:NFMAX*NSUBFMAX)
 *			        Example:
 *			        	For a FACE of type QUAD, there are nine possible values for Indsfh, related to the allowed
 *			        	h-refinements:
 *			        		None       (1 -> 1): Indfh = [0]
 *			        		Isotropic  (1 -> 4): Indfh = [1,4]
 *			        		Horizontal (1 -> 2): Indfh = [5,6]
 *			        		Vertical   (1 -> 2): Indfh = [7,8]
 *
 *			        		Indfh = 0 : Indsfh = 0 (This is a conforming FACE)
 *			        		        1 : Indsfh = 0 (This is the first  sub-face of the h-refined macro FACE)
 *			        		        2 : Indsfh = 1 (This is the second sub-face of the h-refined macro FACE)
 *			        		        3 : Indsfh = 2 (This is the third  sub-face of the h-refined macro FACE)
 *			        		        4 : Indsfh = 3 (This is the fourth sub-face of the h-refined macro FACE)
 *			        		        5 : Indsfh = 0
 *			        		        6 : Indsfh = 1
 *			        		        7 : Indsfh = 0
 *			        		        8 : Indsfh = 1
 *
 *			(Ord) : (Ord)ering. Specifies the index for the ordering between FACEs as seen from adjacent VOLUMEs. When
 *			        interpolating values from VOLUMEs to FACEs, it is implicitly assumed that the FACE is in a reference
 *			        configuration of the current VOLUME. When these interpolated values must be used in relation to the
 *			        neighbouring FACE, their orientation is likely to be incorrect and the reordering necessary to be
 *			        seen correctly is required. This variable provides the index to the ordering array which transfers
 *			        the current ordering to that required as if the FACE were seen from the opposite VOLUME.
 *
 *		VOLUME : Pointer to neighbouring VOLUME.
 *
 *
 */

#include <complex.h>

#include "S_VOLUME.h"

struct S_FACE {
	// Structures
	unsigned int P, type, VfL, VfR, indexg, BC, IndOrdLR, IndOrdRL, level, update, adapt_type;

	// Geometry
	unsigned int curved, typeInt;
	double       *XYZ_fI, *XYZ_fS, *n_fI, *n_fS, *detJF_fI, *detJF_fS, *detJVIn_fI, *detJVOut_fI;

	// Initialization
	unsigned int NvnS;

	// Solving
	char         CDG2_side;
	unsigned int IndA, Boundary;
	double       *RHSL, *RHSR, *LHSLL, *LHSRL, *LHSLR, *LHSRR,
	             **QhatL, **QhatR, **QhatL_WhatL, **QhatL_WhatR, **QhatR_WhatL, **QhatR_WhatR;

	// Linearization testing
	double complex *RHSL_c, *RHSR_c, **QhatL_c, **QhatR_c;

	// structs
	struct S_VOLUME *VL, *VR;
	struct S_FACE  *next, *grpnext, *child0, *parent;

	// ToBeModified (Use in updated versions of the code)
	struct S_SIDE_DATA {
		unsigned int   Indfh, Indmf, Indlfh, Indsfh, IndOrd;
// Likely use S_MULTI_ARRAYs for the variables below (ToBeDeleted)
//		double         **Qhat, **dQhat_dWhatL, **dQhat_dWhatR;
//		double complex **Qhat_c;

		struct S_VOLUME *VOLUME;
	} data[2];
};

extern void update_data_FACE (struct S_FACE *const FACE);

#endif // DPG__S_FACE_h__INCLUDED
