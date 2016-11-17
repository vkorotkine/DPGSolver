// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "element_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "memory_constructors.h"

/*
 *	Purpose:
 *		Provide simple element-related functions: (ToBeModified)
 *			int        is_ELEMENT_present(const unsigned int type);
 *			*S_ELEMENT get_ELEMENT_type(const unsigned int type);
 *			*S_ELEMENT get_ELEMENT_Eclass(const unsigned int Eclass, const unsigned int Esubclass);
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void initialize_ELEMENTs(void)
{
	/*
	 *	Purpose:
	 *		Initialize ELEMENT structures.
	 *
	 *	Comments:
	 *		VeCGmsh is based on converting the gmsh node ordering to the code's ordering convention; for the gmsh
	 *		ordering, see the gmsh manual (Node ordering section). VeE and VeF are then set based on the standard node
	 *		ordering in the code.
	 *
	 *	Notation:
	 *		type     : Flag for element type
	 *		type_h   : Type of ELEMENTs present in h-refined ELEMENT
	 *		Eclass   : (E)lement (class)
	 *		d        : (d)imension
	 *		Nve      : (N)umber of local (ve)rtices
	 *		Nfve     : (N)umber of local (f)ace (ve)rtices
	 *		Nf       : (N)umber of local (f)aces
	 *		NEhref   : (N)umber of (E)lement types after (h)-(ref)inement
	 *		VeCGmsh  : (Ve)rtices of the (C)orners                      (Gmsh ordering)
	 *		VeEcon   : (Ve)rtices of the (E)dges which are (con)forming (Standard ordering)
	 *		VeFcon   : (Ve)rtices of the (F)aces which are (con)forming (Standard ordering)
	 *		Nvref    : (N)umber of (v)olume h-refinements needed for standard refinement.
	 *		NvrefSF  : (N)umber of (v)olume h-refinements needed for (S)um (F)actorized operator application.
	 *		Nfref    : (N)umber of (f)ace h-refinements needed for standard refinement.
	 *
	 *	References:
	 *		http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
	 *
	 */

	// Initialize DB Parameters
	unsigned int Adapt         = DB.Adapt,
	             TETrefineType = DB.TETrefineType;

	// Standard datatypes
	unsigned int type, f, Nf, IndFType;

	struct S_ELEMENT *ELEMENT, *ELEMENT_FACE;

	// POINT
	ELEMENT = New_ELEMENT();
	DB.ELEMENT = ELEMENT;

	ELEMENT->type      = POINT;
	ELEMENT->Eclass    = C_TP;
	ELEMENT->d         = 0;
	ELEMENT->Nve       = 1;
	ELEMENT->Ne        = 1;
	ELEMENT->Nf        = 1;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = POINT;
	ELEMENT->Neve      = 1;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 1;
	ELEMENT->VeCGmsh[0]    = 0;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 0;

	ELEMENT->next = New_ELEMENT();

	// LINE
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = LINE;
	ELEMENT->Eclass    = C_TP;
	ELEMENT->d         = 1;
	ELEMENT->Nve       = 2;
	ELEMENT->Ne        = 2;
	ELEMENT->Nf        = 2;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = LINE;
	ELEMENT->Neve      = 1;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 1; ELEMENT->Nfve[1]   = 1;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1] =    1;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 0;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 1;
	ELEMENT->NrefV[0] = 2;

	ELEMENT->next = New_ELEMENT();

	// TRI
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = TRI;
	ELEMENT->Eclass    = C_SI;
	ELEMENT->d         = 2;
	ELEMENT->Nve       = 3;
	ELEMENT->Ne        = 3;
	ELEMENT->Nf        = 3;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = TRI;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 2; ELEMENT->Nfve[1]   = 2; ELEMENT->Nfve[2]   = 2;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2] =    2;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 1; ELEMENT->VeFcon[0*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 0; ELEMENT->VeFcon[1*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1;
	ELEMENT->NrefV[0] = 4;

	ELEMENT->next = New_ELEMENT();

	// QUAD
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = QUAD;
	ELEMENT->Eclass    = C_TP;
	ELEMENT->d         = 2;
	ELEMENT->Nve       = 4;
	ELEMENT->Ne        = 4;
	ELEMENT->Nf        = 4;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = QUAD;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 2; ELEMENT->Nfve[1]   = 2; ELEMENT->Nfve[2]   = 2; ELEMENT->Nfve[3]   = 2;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1] =    1; ELEMENT->VeCGmsh[2] =    3; ELEMENT->VeCGmsh[3] =    2;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1;
	ELEMENT->VeFcon[3*NFVEMAX  ] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; // ToBeDeleted: Modified from matlab code
	ELEMENT->NrefV[0] = 4; ELEMENT->NrefV[1] = 2; ELEMENT->NrefV[2] = 2;

	ELEMENT->next = New_ELEMENT();

	// TET
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = TET;
	ELEMENT->Eclass    = C_SI;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 4;
	ELEMENT->Ne        = 6;
	ELEMENT->Nf        = 4;
	ELEMENT->NEhref    = 2;
	ELEMENT->type_h[0] = TET; ELEMENT->type_h[1] = PYR;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 3; ELEMENT->Nfve[1]   = 3; ELEMENT->Nfve[2]   = 3; ELEMENT->Nfve[3]   = 3;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 2; ELEMENT->VeCGmsh[3]    = 3;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 2; ELEMENT->VeFcon[0*NFVEMAX+1] = 3; ELEMENT->VeFcon[0*NFVEMAX+2] = 1;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 2; ELEMENT->VeFcon[1*NFVEMAX+1] = 3; ELEMENT->VeFcon[1*NFVEMAX+2] = 0;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 3;
	ELEMENT->VeFcon[3*NFVEMAX  ] = 0; ELEMENT->VeFcon[3*NFVEMAX+1] = 1; ELEMENT->VeFcon[3*NFVEMAX+2] = 2;
	ELEMENT->VeEcon[0*NEVEMAX  ] = 0; ELEMENT->VeEcon[0*NEVEMAX+1] = 1;
	ELEMENT->VeEcon[1*NEVEMAX  ] = 0; ELEMENT->VeEcon[1*NEVEMAX+1] = 2;
	ELEMENT->VeEcon[2*NEVEMAX  ] = 0; ELEMENT->VeEcon[2*NEVEMAX+1] = 3;
	ELEMENT->VeEcon[3*NEVEMAX  ] = 1; ELEMENT->VeEcon[3*NEVEMAX+1] = 2;
	ELEMENT->VeEcon[4*NEVEMAX  ] = 1; ELEMENT->VeEcon[4*NEVEMAX+1] = 3;
	ELEMENT->VeEcon[5*NEVEMAX  ] = 2; ELEMENT->VeEcon[5*NEVEMAX+1] = 3;
	ELEMENT->NrefV[0] = 6;

	ELEMENT->next = New_ELEMENT();

	// HEX
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = HEX;
	ELEMENT->Eclass    = C_TP;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 8;
	ELEMENT->Ne        = 12;
	ELEMENT->Nf        = 6;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = HEX;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 4; ELEMENT->Nfve[1]   = 4; ELEMENT->Nfve[2]   = 4;
	ELEMENT->Nfve[3]   = 4; ELEMENT->Nfve[4]   = 4; ELEMENT->Nfve[5]   = 4;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 3; ELEMENT->VeCGmsh[3]    = 2;
	ELEMENT->VeCGmsh[4]    = 4; ELEMENT->VeCGmsh[5]    = 5; ELEMENT->VeCGmsh[6]    = 7; ELEMENT->VeCGmsh[7]    = 6;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4; ELEMENT->VeFcon[0*NFVEMAX+3] = 6;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3; ELEMENT->VeFcon[1*NFVEMAX+2] = 5; ELEMENT->VeFcon[1*NFVEMAX+3] = 7;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 4; ELEMENT->VeFcon[2*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[3*NFVEMAX  ] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; ELEMENT->VeFcon[3*NFVEMAX+2] = 6; ELEMENT->VeFcon[3*NFVEMAX+3] = 7;
	ELEMENT->VeFcon[4*NFVEMAX  ] = 0; ELEMENT->VeFcon[4*NFVEMAX+1] = 1; ELEMENT->VeFcon[4*NFVEMAX+2] = 2; ELEMENT->VeFcon[4*NFVEMAX+3] = 3;
	ELEMENT->VeFcon[5*NFVEMAX  ] = 4; ELEMENT->VeFcon[5*NFVEMAX+1] = 5; ELEMENT->VeFcon[5*NFVEMAX+2] = 6; ELEMENT->VeFcon[5*NFVEMAX+3] = 7;
	ELEMENT->VeEcon[0 *NEVEMAX  ] = 0; ELEMENT->VeEcon[0 *NEVEMAX+1] = 1;
	ELEMENT->VeEcon[1 *NEVEMAX  ] = 2; ELEMENT->VeEcon[1 *NEVEMAX+1] = 3;
	ELEMENT->VeEcon[2 *NEVEMAX  ] = 4; ELEMENT->VeEcon[2 *NEVEMAX+1] = 5;
	ELEMENT->VeEcon[3 *NEVEMAX  ] = 6; ELEMENT->VeEcon[3 *NEVEMAX+1] = 7;
	ELEMENT->VeEcon[4 *NEVEMAX  ] = 0; ELEMENT->VeEcon[4 *NEVEMAX+1] = 2;
	ELEMENT->VeEcon[5 *NEVEMAX  ] = 1; ELEMENT->VeEcon[5 *NEVEMAX+1] = 3;
	ELEMENT->VeEcon[6 *NEVEMAX  ] = 4; ELEMENT->VeEcon[6 *NEVEMAX+1] = 6;
	ELEMENT->VeEcon[7 *NEVEMAX  ] = 5; ELEMENT->VeEcon[7 *NEVEMAX+1] = 7;
	ELEMENT->VeEcon[8 *NEVEMAX  ] = 0; ELEMENT->VeEcon[8 *NEVEMAX+1] = 4;
	ELEMENT->VeEcon[9 *NEVEMAX  ] = 1; ELEMENT->VeEcon[9 *NEVEMAX+1] = 5;
	ELEMENT->VeEcon[10*NEVEMAX  ] = 2; ELEMENT->VeEcon[10*NEVEMAX+1] = 6;
	ELEMENT->VeEcon[11*NEVEMAX  ] = 3; ELEMENT->VeEcon[11*NEVEMAX+1] = 7;
	ELEMENT->NrefV[0] = 8; ELEMENT->NrefV[1] = 4; ELEMENT->NrefV[2] = 4; ELEMENT->NrefV[3] = 4;
	ELEMENT->NrefV[4] = 2; ELEMENT->NrefV[5] = 2; ELEMENT->NrefV[6] = 2;

	ELEMENT->next = New_ELEMENT();

	// WEDGE
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = WEDGE;
	ELEMENT->Eclass    = C_WEDGE;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 6;
	ELEMENT->Ne        = 9;
	ELEMENT->Nf        = 5;
	ELEMENT->NEhref    = 1;
	ELEMENT->type_h[0] = WEDGE;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 4; ELEMENT->Nfve[1]   = 4; ELEMENT->Nfve[2]   = 4;
	ELEMENT->Nfve[3]   = 3; ELEMENT->Nfve[4]   = 3;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 2;
	ELEMENT->VeCGmsh[3]    = 3; ELEMENT->VeCGmsh[4]    = 4; ELEMENT->VeCGmsh[5]    = 5;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 1; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4; ELEMENT->VeFcon[0*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 0; ELEMENT->VeFcon[1*NFVEMAX+1] = 2; ELEMENT->VeFcon[1*NFVEMAX+2] = 3; ELEMENT->VeFcon[1*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 3; ELEMENT->VeFcon[2*NFVEMAX+3] = 4;
	ELEMENT->VeFcon[3*NFVEMAX  ] = 0; ELEMENT->VeFcon[3*NFVEMAX+1] = 1; ELEMENT->VeFcon[3*NFVEMAX+2] = 2;
	ELEMENT->VeFcon[4*NFVEMAX  ] = 3; ELEMENT->VeFcon[4*NFVEMAX+1] = 4; ELEMENT->VeFcon[4*NFVEMAX+2] = 5;
	ELEMENT->NrefV[0] = 8; ELEMENT->NrefV[1] = 4; ELEMENT->NrefV[2] = 2;

	ELEMENT->next = New_ELEMENT();

	// PYR
	ELEMENT = ELEMENT->next;

	ELEMENT->type      = PYR;
	ELEMENT->Eclass    = C_PYR;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 5;
	ELEMENT->Ne        = 8;
	ELEMENT->Nf        = 5;
	ELEMENT->NEhref    = 2;
	ELEMENT->type_h[0] = PYR; ELEMENT->type_h[1] = TET;
	ELEMENT->Neve      = 2;
	ELEMENT->Neref     = NEREFMAX;
	ELEMENT->Nfve[0]   = 3; ELEMENT->Nfve[1]   = 3; ELEMENT->Nfve[2]   = 3; ELEMENT->Nfve[3]   = 3;
	ELEMENT->Nfve[4]   = 4;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 3; ELEMENT->VeCGmsh[3]   = 2;
	ELEMENT->VeCGmsh[4]    = 4;
	ELEMENT->VeFcon[0*NFVEMAX  ] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[1*NFVEMAX  ] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3; ELEMENT->VeFcon[1*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[2*NFVEMAX  ] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[3*NFVEMAX  ] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; ELEMENT->VeFcon[3*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[4*NFVEMAX  ] = 0; ELEMENT->VeFcon[4*NFVEMAX+1] = 1; ELEMENT->VeFcon[4*NFVEMAX+2] = 2; ELEMENT->VeFcon[4*NFVEMAX+3] = 3;
	ELEMENT->VeEcon[0*NEVEMAX  ] = 2; ELEMENT->VeEcon[0*NEVEMAX+1] = 4;
	ELEMENT->VeEcon[1*NEVEMAX  ] = 0; ELEMENT->VeEcon[1*NEVEMAX+1] = 4;
	ELEMENT->VeEcon[2*NEVEMAX  ] = 0; ELEMENT->VeEcon[2*NEVEMAX+1] = 2;
	ELEMENT->VeEcon[3*NEVEMAX  ] = 3; ELEMENT->VeEcon[3*NEVEMAX+1] = 4;
	ELEMENT->VeEcon[4*NEVEMAX  ] = 1; ELEMENT->VeEcon[4*NEVEMAX+1] = 4;
	ELEMENT->VeEcon[5*NEVEMAX  ] = 1; ELEMENT->VeEcon[5*NEVEMAX+1] = 3;
	ELEMENT->VeEcon[6*NEVEMAX  ] = 0; ELEMENT->VeEcon[6*NEVEMAX+1] = 1;
	ELEMENT->VeEcon[7*NEVEMAX  ] = 2; ELEMENT->VeEcon[7*NEVEMAX+1] = 3;
	ELEMENT->NrefV[0] = 10;

	// No additional ELEMENTs

	// Set pointers for ELEMENT classes and ELEMENT FACEs
	for (ELEMENT = DB.ELEMENT; ELEMENT; ELEMENT = ELEMENT->next) {
		type = ELEMENT->type;
		if (type == POINT) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
		} else if (type == LINE || type == QUAD || type == HEX || type == TRI || type == TET) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACE[0] = get_ELEMENT_FACE(ELEMENT->type,0);
		} else if (type == PYR) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACE[0] = get_ELEMENT_FACE(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACE[1] = get_ELEMENT_FACE(ELEMENT->type,1);
		} else if (type == WEDGE) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENTclass[1]  = get_ELEMENT_Eclass(ELEMENT->type,1);
			ELEMENT->ELEMENT_FACE[0] = get_ELEMENT_FACE(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACE[1] = get_ELEMENT_FACE(ELEMENT->type,1);
		}
	}

	// Set Nvref/Nfref
	switch (Adapt) {
	case ADAPT_0:
	case ADAPT_P:
		for (ELEMENT = DB.ELEMENT; ELEMENT; ELEMENT = ELEMENT->next) {
			ELEMENT->Nvref   = 1;
			ELEMENT->NvrefSF = 1;
		}
		break;
	default: // ADAPT_H and ADAPT_HP
		for (ELEMENT = DB.ELEMENT; ELEMENT; ELEMENT = ELEMENT->next) {
			type = ELEMENT->type;
			switch (type) {
			case POINT:
				ELEMENT->Nvref   = NREFMAXPOINT;
				ELEMENT->NvrefSF = 0; // Not used.
				break;
			case LINE:
				ELEMENT->Nvref   = NREFMAXLINE;
				ELEMENT->NvrefSF = NREFMAXLINE;
				break;
			case TRI:
				ELEMENT->Nvref   = NREFMAXTRI;
				ELEMENT->NvrefSF = NREFMAXTRI;
				break;
			case QUAD:
				ELEMENT->Nvref   = NREFMAXQUAD;
				ELEMENT->NvrefSF = 0; // Not used.
				break;
			case TET:
				if (TETrefineType == TET8)
					ELEMENT->Nvref = 8+1;
				else if (TETrefineType == TET12)
					ELEMENT->Nvref = NREFMAXTET;
				else if (TETrefineType == TET6)
					ELEMENT->Nvref = 6+1;
				else
					printf("Error: Unsupported.\n"), EXIT_MSG;

				ELEMENT->NvrefSF = 1;
				break;
			case HEX:
				ELEMENT->Nvref   = NREFMAXHEX;
				ELEMENT->NvrefSF = 0; // Not used.
				break;
			case WEDGE:
				ELEMENT->Nvref   = NREFMAXWEDGE;
				ELEMENT->NvrefSF = 0; // Not used.
				break;
			case PYR:
				ELEMENT->Nvref   = NREFMAXPYR;
				ELEMENT->NvrefSF = 1;
				break;
			default:
				printf("Error: Unsupported type in Nvref initialization.\n"), exit(1);
				break;
			}
		}
		break;
	}

	for (ELEMENT = DB.ELEMENT->next; ELEMENT; ELEMENT = ELEMENT->next) {
		Nf = ELEMENT->Nf;
		for (f = 0; f < Nf; f++) {
			IndFType = get_IndFType(ELEMENT->Eclass,f);
			ELEMENT_FACE = get_ELEMENT_FACE(ELEMENT->type,IndFType);
			ELEMENT->Nfref[f] = ELEMENT_FACE->Nvref;
		}
	}
}

void finalize_ELEMENTs(void)
{
	/*
	 *	Purpose:
	 *		Finalize ELEMENT structures (initialization).
	 *
	 *	Comments:
	 *		As Gmsh only stores "physical" elements in the mesh file, it is possible that TRI/QUAD ELEMENTs be marked as
	 *		not present when using 3D mixed meshes if there are none on the boundary. Thus, TRI and QUAD ELEMENT
	 *		presence is treated as a special case.
	 *
	 *	Notation:
	 *		present  : Indicator of presence of this element type
	 *		           Options: 0, 1
	 *		NfMax    : (Max)imum (N)umber of (f)aces on an element
	 *		NfveMax  : (Max)imum (N)umber of (f)ace (ve)rtices on an element
	 *		NveMax   : (Max)imum (N)umber of (ve)rtices on an element
	 *		NfrefMax : (Max)imum (N)umber of (f)ace (ref)inements
	 *
	 *	References:
	 */

	// Initialize DB Parameters
	const unsigned int d       = DB.d,
	                   NETotal = DB.NETotal,
	                   *EType  = DB.EType;

	// Standard datatypes
	unsigned int i, NfMax, NfveMax, NveMax, NfrefMax, TRIpresent3D, QUADpresent3D, PYRpresent, type;

	struct S_ELEMENT *ELEMENT;

	TRIpresent3D  = 0;
	QUADpresent3D = 0;
	PYRpresent    = 0;
	for (ELEMENT = DB.ELEMENT; ELEMENT; ELEMENT = ELEMENT->next) {
		type = ELEMENT->type;

		for (i = 0; i < NETotal; i++) {
			if (type == EType[i]) {
				ELEMENT->present = 1;

				if (!TRIpresent3D && (type == TET || type == WEDGE || type == PYR))
					TRIpresent3D = 1;
				if (!QUADpresent3D && (type == HEX || type == WEDGE || type == PYR))
					QUADpresent3D = 1;

				break;
			}
		}
		if (type == PYR && ELEMENT->present == 1)
			PYRpresent = 1;
	}

	// Ensure that 2D ELEMENTs are marked correctly on 3D mixed meshes
	if (!DB.MPIrank && !DB.Testing)
		printf("      ELEMENT types present: ");

	for (ELEMENT = DB.ELEMENT; ELEMENT; ELEMENT = ELEMENT->next) {
		type = ELEMENT->type;
		if ((type == TRI && TRIpresent3D) || (type == QUAD && QUADpresent3D) || (type == TET && PYRpresent))
			ELEMENT->present = 1;

		if (ELEMENT->present) {
			if (!DB.MPIrank && !DB.Testing)
				printf("%d, ",type);
		}
	}
	if (!DB.MPIrank && !DB.Testing)
		printf("\n");

	if (d == 1) {
		NfMax    = 2;
		NfveMax  = 1;
		NveMax   = 2;
	} else if (d == 2) {
		NfMax    = 3;
		NfveMax  = 2;
		NveMax   = 3;

		ELEMENT = get_ELEMENT_type(QUAD);
		if (ELEMENT->present == 1) {
			NfMax  = 4;
			NveMax = 4;
		}
	} else {
		NfMax    = 4;
		NfveMax  = 3;
		NveMax   = 4;

		ELEMENT = get_ELEMENT_type(PYR);
		if (ELEMENT->present == 1) {
			NfMax    = 5;
			NfveMax  = 4;
			NveMax   = 5;
		}

		ELEMENT = get_ELEMENT_type(WEDGE);
		if (ELEMENT->present == 1) {
			NfMax    = 5;
			NfveMax  = 4;
			NveMax   = 6;
		}

		ELEMENT = get_ELEMENT_type(HEX);
		if (ELEMENT->present == 1) {
			NfMax    = 6;
			NfveMax  = 4;
			NveMax   = 8;
		}
	}
	NfrefMax = NFREFMAX;

	DB.NfMax    = NfMax;
	DB.NfveMax  = NfveMax;
	DB.NveMax   = NveMax;
	DB.NfrefMax = NfrefMax;
}

unsigned int is_ELEMENT_present(const unsigned int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT) {
		if (type == ELEMENT->type) {
			if (ELEMENT->present)
				return 1;
			else
				return 0;
		}
		ELEMENT = ELEMENT->next;
	}
	printf("Error: Element type not found (present).\n"), exit(1);
}

unsigned int get_Eclass(const unsigned int type)
{
	if (type == POINT || type == LINE || type == QUAD || type == HEX)
		return C_TP;
	else if (type == TRI || type == TET)
		return C_SI;
	else if (type == PYR)
		return C_PYR;
	else if (type == WEDGE)
		return C_WEDGE;

	printf("Error: There is not yet an element class associated with the type provided.\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_type(const unsigned int type)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	while (ELEMENT) {
		if (type == ELEMENT->type)
			return ELEMENT;

		ELEMENT = ELEMENT->next;
	}
	printf("Error: Element type not found (type).\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_F_type(const unsigned int type, const unsigned int f)
{
	unsigned int IndFType;

	struct S_ELEMENT *ELEMENT_F;

	IndFType = get_IndFType(get_ELEMENT_type(type)->Eclass,f);

	if (type == LINE) {
		ELEMENT_F = get_ELEMENT_type(POINT);
	} else if (type == TRI || type == QUAD) {
		ELEMENT_F = get_ELEMENT_type(LINE);
	} else if (type == TET || (type == WEDGE && IndFType == 1) || (type == PYR && IndFType == 0)) {
		ELEMENT_F = get_ELEMENT_type(TRI);
	} else if (type == HEX || (type == WEDGE && IndFType == 0) || (type == PYR && IndFType == 1)) {
		ELEMENT_F = get_ELEMENT_type(QUAD);
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
	return ELEMENT_F;
}

struct S_ELEMENT *get_ELEMENT_Eclass(const unsigned int type, const unsigned int IndEclass)
{
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	if (type == POINT || type == LINE || type == QUAD || type == HEX || (type == WEDGE && IndEclass == 1)) {
		while (ELEMENT) {
			if (ELEMENT->type == LINE)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (type == TRI || type == TET || (type == WEDGE && IndEclass == 0)) {
		while (ELEMENT) {
			if (ELEMENT->type == TRI)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	} else if (type == PYR) {
		while (ELEMENT) {
			if (ELEMENT->type == PYR)
				return ELEMENT;

			ELEMENT = ELEMENT->next;
		}
	}
	printf("Error: Element class not found.\n"), exit(1);
}

struct S_ELEMENT *get_ELEMENT_FACE(const unsigned int type, const unsigned int IndEclass)
{
	// Quite similar to get_ELEMENT_F_type ... Likely delete one of this function (ToBeDeleted)
	struct S_ELEMENT *ELEMENT = DB.ELEMENT;

	if (type == LINE) {
		while (ELEMENT) {
			if (ELEMENT->type == POINT)
				return ELEMENT;
			ELEMENT = ELEMENT->next;
		}
	} else if (type == TRI || type == QUAD) {
		while (ELEMENT) {
			if (ELEMENT->type == LINE)
				return ELEMENT;
			ELEMENT = ELEMENT->next;
		}
	} else if (type == TET || (type == WEDGE && IndEclass == 1) || (type == PYR && IndEclass == 0)) {
		while (ELEMENT) {
			if (ELEMENT->type == TRI)
				return ELEMENT;
			ELEMENT = ELEMENT->next;
		}
	} else if (type == HEX || (type == WEDGE && IndEclass == 0) || (type == PYR && IndEclass == 1)) {
		while (ELEMENT) {
			if (ELEMENT->type == QUAD)
				return ELEMENT;
			ELEMENT = ELEMENT->next;
		}
	}
	printf("Error: Element FACE of type %d and IndFType %d was not found.\n",type,IndEclass), exit(1);
}

unsigned int get_IndFType(const unsigned int Eclass, const unsigned int f)
{
	switch (Eclass) {
	case C_TP:
	case C_SI:
		return 0;
		break;
	case C_PYR:
		if (f < 4) return 0;
		else       return 1;
		break;
	case C_WEDGE:
		if (f < 3) return 0;
		else       return 1;
		break;
	default:
		printf("Error: Unsupported Eclass/f combination in get_IndFType.\n"), exit(1);
		break;
	}
}
