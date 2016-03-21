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
 *		This function is the analogue of the following original function(s):
 *
 *		Note: 0-based indexing.
 *		Note: The VeC, VeE, VeF ordering is based on the gmsh node return order convention.
 *		      1D: VeC == VeE == VeF
 *		      2D: VeC != VeE == VeF
 *		      3D: VeC != VeE != VeF
 *		          VeE currently unused => Check values later (ToBeDeleted)
 *
 *	Notation:
 *		present : Indicator of presence of this element type
 *		          Options: 0, 1
 *		type    : Flag for element type
 *		d       : (d)imension
 *		Nve     : (N)umber of local (ve)rtices
 *		Nfve    : (N)umber of local (f)acet (ve)rtices
 *		Nf      : (N)umber of local (f)acets
 *		VeC     : (Ve)rtices of the (C)orners
 *		VeE     : (Ve)rtices of the (E)dges
 *		VeF     : (Ve)rtices of the (F)acets
 *
 *	References:
 *
 */

void setup_mesh()
{
	int type, i, NfMax, NfveMax;

	struct S_ELEMENT *ELEMENT;

	// POINT
	ELEMENT = New_ELEMENT();
	DB.ELEMENT = ELEMENT;

	ELEMENT->present = 0;
	ELEMENT->type    = POINT;
	ELEMENT->d       = 0;
	ELEMENT->Nve     = 1;
	ELEMENT->Nf      = 1;
	ELEMENT->Nfve[0] = 1;
	ELEMENT->VeC[0]     = 0;
	ELEMENT->VeE[0*2+0] = 0;
	ELEMENT->VeF[0*4+0] = 0;

	ELEMENT->next = New_ELEMENT();

	// LINE
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = LINE;
	ELEMENT->d       = 1;
	ELEMENT->Nve     = 2;
	ELEMENT->Nf      = 2;
	ELEMENT->Nfve[0] = 1;
	ELEMENT->VeC[0]     = 0; ELEMENT->VeC[1]     = 1;
	ELEMENT->VeE[0*2+0] = 0;
	ELEMENT->VeE[1*2+0] = 1;
	ELEMENT->VeF[0*4+0] = 0;
	ELEMENT->VeF[1*4+0] = 1;

	ELEMENT->next = New_ELEMENT();

	// TRI
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = TRI;
	ELEMENT->d       = 2;
	ELEMENT->Nve     = 3;
	ELEMENT->Nf      = 3;
	ELEMENT->Nfve[0] = 2;
	ELEMENT->VeC[0]     = 1; ELEMENT->VeC[1]     = 2; ELEMENT->VeC[2]     = 0;
	ELEMENT->VeE[0*2+0] = 1; ELEMENT->VeE[0*2+1] = 0;
	ELEMENT->VeE[1*2+0] = 2; ELEMENT->VeE[1*2+1] = 0;
	ELEMENT->VeE[2*2+0] = 1; ELEMENT->VeE[2*2+1] = 2;
	ELEMENT->VeF[0*4+0] = 1; ELEMENT->VeF[0*4+1] = 0;
	ELEMENT->VeF[1*4+0] = 2; ELEMENT->VeF[1*4+1] = 0;
	ELEMENT->VeF[2*4+0] = 1; ELEMENT->VeF[2*4+1] = 2;

	ELEMENT->next = New_ELEMENT();

	// QUAD
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = QUAD;
	ELEMENT->d       = 2;
	ELEMENT->Nve     = 4;
	ELEMENT->Nf      = 4;
	ELEMENT->Nfve[0] = 2;
	ELEMENT->VeC[0]     = 0; ELEMENT->VeC[1]     = 1; ELEMENT->VeC[2]     = 3; ELEMENT->VeC[3]     = 2;
	ELEMENT->VeE[0*2+0] = 0; ELEMENT->VeE[0*2+1] = 3;
	ELEMENT->VeE[1*2+0] = 1; ELEMENT->VeE[1*2+1] = 2;
	ELEMENT->VeE[2*2+0] = 0; ELEMENT->VeE[2*2+1] = 1;
	ELEMENT->VeE[3*2+0] = 2; ELEMENT->VeE[3*2+1] = 3;
	ELEMENT->VeF[0*4+0] = 0; ELEMENT->VeF[0*4+1] = 3;
	ELEMENT->VeF[1*4+0] = 1; ELEMENT->VeF[1*4+1] = 2;
	ELEMENT->VeF[2*4+0] = 0; ELEMENT->VeF[2*4+1] = 1;
	ELEMENT->VeF[3*4+0] = 2; ELEMENT->VeF[3*4+1] = 3;

	ELEMENT->next = New_ELEMENT();

	// TET
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = TET;
	ELEMENT->d       = 3;
	ELEMENT->Nve     = 4;
	ELEMENT->Nf      = 4;
	ELEMENT->Nfve[0] = 3;
	ELEMENT->VeC[0]     = 1; ELEMENT->VeC[1]     = 3; ELEMENT->VeC[2]     = 2; ELEMENT->VeC[3]     = 0;
	ELEMENT->VeE[0*2+0] = 0; ELEMENT->VeE[0*2+1] = 2;
	ELEMENT->VeE[1*2+0] = 2; ELEMENT->VeE[1*2+1] = 1;
	ELEMENT->VeE[2*2+0] = 1; ELEMENT->VeE[2*2+1] = 0;
	ELEMENT->VeE[3*2+0] = 0; ELEMENT->VeE[3*2+1] = 3;
	ELEMENT->VeE[4*2+0] = 3; ELEMENT->VeE[4*2+1] = 2;
	ELEMENT->VeE[5*2+0] = 1; ELEMENT->VeE[5*2+1] = 3;
	ELEMENT->VeF[0*4+0] = 0; ELEMENT->VeF[0*4+1] = 2; ELEMENT->VeF[0*4+2] = 1;
	ELEMENT->VeF[1*4+0] = 0; ELEMENT->VeF[1*4+1] = 3; ELEMENT->VeF[1*4+2] = 2;
	ELEMENT->VeF[2*4+0] = 0; ELEMENT->VeF[2*4+1] = 1; ELEMENT->VeF[2*4+2] = 3;
	ELEMENT->VeF[3*4+0] = 2; ELEMENT->VeF[3*4+1] = 1; ELEMENT->VeF[3*4+2] = 3;

	ELEMENT->next = New_ELEMENT();

	// HEX
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = HEX;
	ELEMENT->d       = 3;
	ELEMENT->Nve     = 8;
	ELEMENT->Nf      = 6;
	ELEMENT->Nfve[0] = 4;
	ELEMENT->VeC[0]      = 0; ELEMENT->VeC[1]      = 1; ELEMENT->VeC[2]      = 3; ELEMENT->VeC[3]      = 2;
	ELEMENT->VeC[4]      = 4; ELEMENT->VeC[5]      = 5; ELEMENT->VeC[6]      = 7; ELEMENT->VeC[7]      = 6;
	ELEMENT->VeE[0*2+0]  = 0; ELEMENT->VeE[0*2+1]  = 1;
	ELEMENT->VeE[1*2+0]  = 1; ELEMENT->VeE[1*2+1]  = 2;
	ELEMENT->VeE[2*2+0]  = 2; ELEMENT->VeE[2*2+1]  = 3;
	ELEMENT->VeE[3*2+0]  = 3; ELEMENT->VeE[3*2+1]  = 0;
	ELEMENT->VeE[4*2+0]  = 4; ELEMENT->VeE[4*2+1]  = 5;
	ELEMENT->VeE[5*2+0]  = 5; ELEMENT->VeE[5*2+1]  = 6;
	ELEMENT->VeE[6*2+0]  = 6; ELEMENT->VeE[6*2+1]  = 7;
	ELEMENT->VeE[7*2+0]  = 7; ELEMENT->VeE[7*2+1]  = 4;
	ELEMENT->VeE[8*2+0]  = 0; ELEMENT->VeE[8*2+1]  = 4;
	ELEMENT->VeE[9*2+0]  = 1; ELEMENT->VeE[9*2+1]  = 5;
	ELEMENT->VeE[10*2+0] = 2; ELEMENT->VeE[10*2+1] = 6;
	ELEMENT->VeE[11*2+0] = 3; ELEMENT->VeE[11*2+1] = 7;
	ELEMENT->VeF[0*4+0]  = 0; ELEMENT->VeF[0*4+1]  = 4; ELEMENT->VeF[0*4+2]  = 7; ELEMENT->VeF[0*4+3]  = 3;
	ELEMENT->VeF[1*4+0]  = 5; ELEMENT->VeF[1*4+1]  = 1; ELEMENT->VeF[1*4+2]  = 2; ELEMENT->VeF[1*4+3]  = 6;
	ELEMENT->VeF[2*4+0]  = 1; ELEMENT->VeF[2*4+1]  = 5; ELEMENT->VeF[2*4+2]  = 4; ELEMENT->VeF[2*4+3]  = 0;
	ELEMENT->VeF[3*4+0]  = 6; ELEMENT->VeF[3*4+1]  = 2; ELEMENT->VeF[3*4+2]  = 3; ELEMENT->VeF[3*4+3]  = 7;
	ELEMENT->VeF[4*4+0]  = 2; ELEMENT->VeF[4*4+1]  = 1; ELEMENT->VeF[4*4+2]  = 0; ELEMENT->VeF[4*4+3]  = 3;
	ELEMENT->VeF[5*4+0]  = 5; ELEMENT->VeF[5*4+1]  = 6; ELEMENT->VeF[5*4+2]  = 7; ELEMENT->VeF[5*4+3]  = 4;

	ELEMENT->next = New_ELEMENT();

	// WEDGE (ToBeModified)
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = WEDGE;
	ELEMENT->d       = 3;
	ELEMENT->Nve     = 6;
	ELEMENT->Nf      = 5;
	ELEMENT->Nfve[0] = 3; ELEMENT->Nfve[1] = 4;

	ELEMENT->next = New_ELEMENT();

	// PYR (ToBeModified)
	ELEMENT = ELEMENT->next;

	ELEMENT->present = 0;
	ELEMENT->type    = PYR;
	ELEMENT->d       = 3;
	ELEMENT->Nve     = 5;
	ELEMENT->Nf      = 5;
	ELEMENT->Nfve[0] = 3; ELEMENT->Nfve[1] = 4;

	// No additional ELEMENTs

	// Read mesh file
	if (!DB.MPIrank)
		printf("    Read MeshFile\n");
	gmsh_reader();

	// Initialize DB Parameters set in gmsh_reader.c
	int NETotal = DB.NETotal,
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
		NfMax = 2;
		NfveMax = 1;
	} else if (DB.d == 2) {
		NfMax = 3;
		NfveMax = 2;

		ELEMENT = DB.ELEMENT; while (ELEMENT->type != QUAD) ELEMENT = ELEMENT->next;
		if (ELEMENT->present == 1)
			NfMax = 4;
	} else {
		NfMax = 4;
		NfveMax = 3;

		ELEMENT = DB.ELEMENT; while (ELEMENT->type != WEDGE) ELEMENT = ELEMENT->next;
		if (ELEMENT->present == 1) {
			NfMax = 5;
			NfveMax = 4;
		}

		ELEMENT = DB.ELEMENT; while (ELEMENT->type != PYR) ELEMENT = ELEMENT->next;
		if (ELEMENT->present == 1) {
			NfMax = 5;
			NfveMax = 4;
		}

		ELEMENT = DB.ELEMENT; while (ELEMENT->type != HEX) ELEMENT = ELEMENT->next;
		if (ELEMENT->present == 1) {
			NfMax = 6;
			NfveMax = 4;
		}
	}

	DB.NfMax   = NfMax;
	DB.NfveMax = NfveMax;

	// Max Order Check
	for (ELEMENT = DB.ELEMENT; ELEMENT != NULL; ELEMENT = ELEMENT->next) {
		if (ELEMENT->type == TRI && ELEMENT->present == 1) {
			if (DB.PFrc[DB.NP][1] > 8)
				printf("\n      WARNING: Reduce PFr if using 2D WS nodes (Max: P8).\n\n");
		} else if (ELEMENT->type == TET && ELEMENT->present == 1) {
			if (DB.PFrc[DB.NP][1] > 6)
				printf("\n      WARNING: Reduce PFr if using 3D WS nodes (Max: P6).\n\n");
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
