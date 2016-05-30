// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up mesh related parameters and initialize the ELEMENT structures.
 *
 *	Comments:
 *		0-based indexing.
 *		VeCGmsh is based on converting the gmsh node ordering to the code's ordering convention; for the gmsh ordering,
 *		see	the gmsh manual (Node ordering section). VeE and VeF are then set based on the standard node ordering in
 *		the code.
 *
 *	Notation:
 *		present  : Indicator of presence of this element type
 *		           Options: 0, 1
 *		type     : Flag for element type
 *		d        : (d)imension
 *		Nve      : (N)umber of local (ve)rtices
 *		Nfve     : (N)umber of local (f)acet (ve)rtices
 *		Nf       : (N)umber of local (f)acets
 *		VeCGmsh  : (Ve)rtices of the (C)orners                       (Gmsh ordering)
 *		VeFcon   : (Ve)rtices of the (F)acets which are (con)forming (Standard ordering)
 *
 *		NfMax    : (Max)imum (N)umber of (f)acets on an element
 *		NfveMax  : (Max)imum (N)umber of (f)acet (ve)rtices on an element
 *		NveMax   : (Max)imum (N)umber of (ve)rtices on an element
 *		NfrefMax : (Max)imum (N)umber of (f)acet (ref)inements
 *
 *	References:
 *		http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
 */

void setup_mesh()
{
	unsigned int i, type, NfveMax, NfMax, NveMax, NfrefMax;

	struct S_ELEMENT *ELEMENT;

	// POINT
	ELEMENT = New_ELEMENT();
	DB.ELEMENT = ELEMENT;

	ELEMENT->present   = 0;
	ELEMENT->type      = POINT;
	ELEMENT->d         = 0;
	ELEMENT->Nve       = 1;
	ELEMENT->Nf        = 1;
	ELEMENT->Nfve[0]   = 1;
	ELEMENT->VeCGmsh[0]    = 0;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 0;

	ELEMENT->next = New_ELEMENT();

	// LINE
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = LINE;
	ELEMENT->d         = 1;
	ELEMENT->Nve       = 2;
	ELEMENT->Nf        = 2;
	ELEMENT->Nfve[0]   = 1; ELEMENT->Nfve[1]   = 1;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1] =    1;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 0;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 1;

	ELEMENT->next = New_ELEMENT();

	// TRI
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = TRI;
	ELEMENT->d         = 2;
	ELEMENT->Nve       = 3;
	ELEMENT->Nf        = 3;
	ELEMENT->Nfve[0]   = 2; ELEMENT->Nfve[1]   = 2; ELEMENT->Nfve[2]   = 2;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2] =    2;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 1; ELEMENT->VeFcon[0*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 0; ELEMENT->VeFcon[1*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1;

	ELEMENT->next = New_ELEMENT();

	// QUAD
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = QUAD;
	ELEMENT->d         = 2;
	ELEMENT->Nve       = 4;
	ELEMENT->Nf        = 4;
	ELEMENT->Nfve[0]   = 2; ELEMENT->Nfve[1]   = 2; ELEMENT->Nfve[2]   = 2; ELEMENT->Nfve[3]   = 2;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1] =    1; ELEMENT->VeCGmsh[2] =    3; ELEMENT->VeCGmsh[3] =    2;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1;
	ELEMENT->VeFcon[3*NFVEMAX+0] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; // ToBeDeleted: Modified from matlab code

	ELEMENT->next = New_ELEMENT();

	// TET
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = TET;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 4;
	ELEMENT->Nf        = 4;
	ELEMENT->Nfve[0]   = 3; ELEMENT->Nfve[1]   = 3; ELEMENT->Nfve[2]   = 3; ELEMENT->Nfve[3]   = 3;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 2; ELEMENT->VeCGmsh[3]    = 3;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 1; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 3;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 0; ELEMENT->VeFcon[1*NFVEMAX+1] = 2; ELEMENT->VeFcon[1*NFVEMAX+2] = 3;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 3;
	ELEMENT->VeFcon[3*NFVEMAX+0] = 0; ELEMENT->VeFcon[3*NFVEMAX+1] = 1; ELEMENT->VeFcon[3*NFVEMAX+2] = 2;

	ELEMENT->next = New_ELEMENT();

	// HEX
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = HEX;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 8;
	ELEMENT->Nf        = 6;
	ELEMENT->Nfve[0]   = 4; ELEMENT->Nfve[1]   = 4; ELEMENT->Nfve[2]   = 4;
	ELEMENT->Nfve[3]   = 4; ELEMENT->Nfve[4]   = 4; ELEMENT->Nfve[5]   = 4;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 3; ELEMENT->VeCGmsh[3]    = 2;
	ELEMENT->VeCGmsh[4]    = 4; ELEMENT->VeCGmsh[5]    = 5; ELEMENT->VeCGmsh[6]    = 7; ELEMENT->VeCGmsh[7]    = 6;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4; ELEMENT->VeFcon[0*NFVEMAX+3] = 6;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3; ELEMENT->VeFcon[1*NFVEMAX+2] = 5; ELEMENT->VeFcon[1*NFVEMAX+3] = 7;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 4; ELEMENT->VeFcon[2*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[3*NFVEMAX+0] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; ELEMENT->VeFcon[3*NFVEMAX+2] = 6; ELEMENT->VeFcon[3*NFVEMAX+3] = 7;
	ELEMENT->VeFcon[4*NFVEMAX+0] = 0; ELEMENT->VeFcon[4*NFVEMAX+1] = 1; ELEMENT->VeFcon[4*NFVEMAX+2] = 2; ELEMENT->VeFcon[4*NFVEMAX+3] = 3;
	ELEMENT->VeFcon[5*NFVEMAX+0] = 4; ELEMENT->VeFcon[5*NFVEMAX+1] = 5; ELEMENT->VeFcon[5*NFVEMAX+2] = 6; ELEMENT->VeFcon[5*NFVEMAX+3] = 7;

	ELEMENT->next = New_ELEMENT();

	// WEDGE (ToBeModified)
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = WEDGE;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 6;
	ELEMENT->Nf        = 5;
	ELEMENT->Nfve[0]   = 4; ELEMENT->Nfve[1]   = 4; ELEMENT->Nfve[2]   = 4;
	ELEMENT->Nfve[3]   = 3; ELEMENT->Nfve[4]   = 3;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 2;
	ELEMENT->VeCGmsh[3]    = 3; ELEMENT->VeCGmsh[4]    = 4; ELEMENT->VeCGmsh[5]    = 5;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 1; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4; ELEMENT->VeFcon[0*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 0; ELEMENT->VeFcon[1*NFVEMAX+1] = 2; ELEMENT->VeFcon[1*NFVEMAX+2] = 3; ELEMENT->VeFcon[1*NFVEMAX+3] = 5;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 3; ELEMENT->VeFcon[2*NFVEMAX+3] = 4;
	ELEMENT->VeFcon[3*NFVEMAX+0] = 0; ELEMENT->VeFcon[3*NFVEMAX+1] = 1; ELEMENT->VeFcon[3*NFVEMAX+2] = 2;
	ELEMENT->VeFcon[4*NFVEMAX+0] = 3; ELEMENT->VeFcon[4*NFVEMAX+1] = 4; ELEMENT->VeFcon[4*NFVEMAX+2] = 5;

	ELEMENT->next = New_ELEMENT();

	// PYR (ToBeModified)
	ELEMENT = ELEMENT->next;

	ELEMENT->present   = 0;
	ELEMENT->type      = PYR;
	ELEMENT->d         = 3;
	ELEMENT->Nve       = 5;
	ELEMENT->Nf        = 5;
	ELEMENT->Nfve[0]   = 3; ELEMENT->Nfve[1]   = 3; ELEMENT->Nfve[2]   = 3; ELEMENT->Nfve[3]   = 3;
	ELEMENT->Nfve[4]   = 4;
	ELEMENT->VeCGmsh[0]    = 0; ELEMENT->VeCGmsh[1]    = 1; ELEMENT->VeCGmsh[2]    = 3; ELEMENT->VeCGmsh[3]   = 2;
	ELEMENT->VeCGmsh[4]    = 4;
	ELEMENT->VeFcon[0*NFVEMAX+0] = 0; ELEMENT->VeFcon[0*NFVEMAX+1] = 2; ELEMENT->VeFcon[0*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[1*NFVEMAX+0] = 1; ELEMENT->VeFcon[1*NFVEMAX+1] = 3; ELEMENT->VeFcon[1*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[2*NFVEMAX+0] = 0; ELEMENT->VeFcon[2*NFVEMAX+1] = 1; ELEMENT->VeFcon[2*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[3*NFVEMAX+0] = 2; ELEMENT->VeFcon[3*NFVEMAX+1] = 3; ELEMENT->VeFcon[3*NFVEMAX+2] = 4;
	ELEMENT->VeFcon[4*NFVEMAX+0] = 0; ELEMENT->VeFcon[4*NFVEMAX+1] = 1; ELEMENT->VeFcon[4*NFVEMAX+2] = 2; ELEMENT->VeFcon[4*NFVEMAX+3] = 3;

	// No additional ELEMENTs

	// Set pointers for ELEMENT classes and ELEMENT FACETs
	for (ELEMENT = DB.ELEMENT; ELEMENT != NULL; ELEMENT = ELEMENT->next) {
		type = ELEMENT->type;
		if (type == POINT) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
		} else if (type == LINE || type == QUAD || type == HEX || type == TRI || type == TET) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACET[0] = get_ELEMENT_FACET(ELEMENT->type,0);
		} else if (type == PYR) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACET[0] = get_ELEMENT_FACET(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACET[1] = get_ELEMENT_FACET(ELEMENT->type,1);
		} else if (type == WEDGE) {
			ELEMENT->ELEMENTclass[0]  = get_ELEMENT_Eclass(ELEMENT->type,0);
			ELEMENT->ELEMENTclass[1]  = get_ELEMENT_Eclass(ELEMENT->type,1);
			ELEMENT->ELEMENT_FACET[0] = get_ELEMENT_FACET(ELEMENT->type,0);
			ELEMENT->ELEMENT_FACET[1] = get_ELEMENT_FACET(ELEMENT->type,1);
		}
	}

	// Read mesh file
	if (!DB.MPIrank)
		printf("    Read MeshFile\n");
	gmsh_reader();

	// Initialize DB Parameters set in gmsh_reader.c
	unsigned int NETotal = DB.NETotal,
	             *EType  = DB.EType;

	for (ELEMENT = DB.ELEMENT; ELEMENT != NULL; ELEMENT = ELEMENT->next) {
		type = ELEMENT->type;
		for (i = 0; i < NETotal; i++) {
			if (type == EType[i]) {
				ELEMENT->present = 1;
				break;
			}
		}
	}

	if (DB.d == 1) {
		NfMax    = 2;
		NfveMax  = 1;
		NveMax   = 2;
	} else if (DB.d == 2) {
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

	// Max Order Check
	for (ELEMENT = DB.ELEMENT; ELEMENT != NULL; ELEMENT = ELEMENT->next) {
		if (ELEMENT->type == TRI && ELEMENT->present == 1) {
			if (DB.PFrc[DB.PMax][1] > 8)
				printf("\n      WARNING: Reduce PFr if using 2D WSH nodes (Max: P8).\n\n");
		} else if (ELEMENT->type == TET && ELEMENT->present == 1) {
			if (DB.PFrc[DB.PMax][1] > 6)
				printf("\n      WARNING: Reduce PFr if using 3D WSH nodes (Max: P6).\n\n");
		}
	}

	// Build Connectivity
	if (!DB.MPIrank)
		printf("    Set up connectivity\n");
	setup_connectivity();

	// Modify connectivity if periodic
	if (DB.NPVe != 0) {
		if (!DB.MPIrank)
			printf("    Modify connectivity for periodic\n");
		setup_periodic();
	}
}
