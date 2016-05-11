#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up VOLUME and FACET structures.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

void setup_structures(void)
{
	// Initialize DB Parameters
	unsigned int d       = DB.d,
	             AC      = DB.AC,
	             P       = DB.P,
	             NP      = DB.NP,
	             NV      = DB.NV,
	             NVC     = DB.NVC,
	             *NE     = DB.NE,
	             *EToVe  = DB.EToVe,
	             *EType  = DB.EType,
	             *EToPrt = DB.EToPrt,
	             *VC     = DB.VC;
	int          MPIrank = DB.MPIrank;
	double       *VeXYZ  = DB.VeXYZ;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int i, iMax, v, dim, ve,
	             IndE, IndVC, IndVgrp,
	             Vs, vlocal, NVlocal, NECgrp, NVgrp,
				 indexg, NvnGs,
	             uMPIrank;
	double       *XYZc;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME *VOLUME, **Vgrp, **Vgrp_tmp;

	// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
	NECgrp = 0;

	uMPIrank = MPIrank;

	if      (d == 1) NECgrp = 1;
	else if (d == 2) NECgrp = 2;
	else if (d == 3) NECgrp = 4;

	NVgrp = NECgrp*NP*2;
	Vgrp     = malloc(NVgrp * sizeof *Vgrp);     // keep
	Vgrp_tmp = malloc(NVgrp * sizeof *Vgrp_tmp); // free
	for (i = 0, iMax = NVgrp; iMax--; i++) {
		Vgrp[i]     = NULL;
		Vgrp_tmp[i] = NULL;
	}

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;

	// Note: only initialize volumes on the current processor.
	for (v = 0, NVlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank)
			NVlocal++;
	}

	for (v = 0, IndE = Vs, IndVC = 0, vlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank) {
			// General
			VOLUME->indexl = vlocal;
			VOLUME->indexg = v;
			VOLUME->P      = P;
			VOLUME->type   = EType[IndE];
			VOLUME->Eclass = get_Eclass(VOLUME->type);

			if (AC || v == VC[IndVC]) {
				VOLUME->curved = 1;
				IndVC++;
			} else {
				VOLUME->curved = 0;
			}

			// Connectivity
			/*	Element connectivity numbering will be based on numbering supporting h-refinement; for this reason it may
			 *	seem overly complex while storing the connectivity information on the initially conforming mesh.
			 *	Supported h-refinements include:
			 *		TET: Isotropic refinement into 4 TET, 2 PYR (PYR orientation tbd)
			 *		HEX: Isotropic and anisotropic refinement of opposite facets into 2 horizontal, 2 vertical, 4 equal
			 *		     QUADs.
			 *		WEDGE: Isotropic refinement of opposing TRI facets, isotropic and anisotropic refinement of all QUAD
			 *		       facets.
			 *		PYR: Isotropic refinement into 6 PYR, 4 TET.
			 *	Max number of faces (including all h-refined variations:
			 *		TET: 4*(1+4) = 20
			 *		HEX: 6*(1+2*2+4) = 54
			 *		WEDGE: 2*(1+4)+3*(1+2*2+4) = 37
			 *		PYR: 4*(1+4)+1*(1+2*2+4) = 25
			 *
			 *			MAX FACES: 6*9 = 54
			 */


			// Geometry
			ELEMENT = get_ELEMENT_type(VOLUME->type);

			if (VOLUME->Eclass == C_TP)
				NvnGs = pow(ELEMENT->ELEMENTclass[0]->NvnGs[0],d);
			else if (VOLUME->Eclass == C_WEDGE)
				NvnGs = (ELEMENT->ELEMENTclass[0]->NvnGs[0])*(ELEMENT->ELEMENTclass[1]->NvnGs[0]);
			else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR)
				NvnGs = ELEMENT->NvnGs[0];
			else
				printf("Error: Unsupported element type setup_struct (NvnGs).\n"), exit(1);

			XYZc = malloc(NvnGs*d * sizeof *XYZc); // keep
			VOLUME->XYZc = XYZc;

			// XYZc may be interpreted as [X Y Z] where each of X, Y, Z are column vectors (ToBeDeleted)
			// Add this comment to notation section.
			indexg = VOLUME->indexg;
			for (ve = 0; ve < NvnGs; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZc[dim*NvnGs+ve] = VeXYZ[EToVe[(Vs+indexg)*8+ve]*d+dim];
			}}

			// MPI
			IndVgrp = ((VOLUME->Eclass)*NP*2)+(VOLUME->P*2)+(VOLUME->curved);
			if (Vgrp[IndVgrp] == NULL)
				Vgrp[IndVgrp] = VOLUME;
			else
				Vgrp_tmp[IndVgrp]->grpnext = VOLUME;

			Vgrp_tmp[IndVgrp] = VOLUME;

			if (vlocal != NVlocal-1) {
				VOLUME->next = New_VOLUME();
				VOLUME = VOLUME->next;
			}

			vlocal++;
		} else {
			// Ensure that appropriate global indices are incremented if necessary
			if (v == VC[IndVC])
				IndVC++;
		}

		IndE++;
	}
	free(Vgrp_tmp);

	if (!AC && IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

/*
for (i = 0, iMax = NVgrp; iMax--; i++) {
	for (VOLUME = Vgrp[i]; VOLUME != NULL; VOLUME = VOLUME->grpnext) {
		printf("%d %d %d %d\n",i,VOLUME->Eclass,VOLUME->P,VOLUME->curved);
	}
	printf("\t\t%p\n",Vgrp[i]);
}
*/
//exit(1);


	// Assign/Overwrite DB parameters
	DB.NV     = NVlocal;
	DB.NECgrp = NECgrp;

	DB.Vgrp = Vgrp;

//VOLUME = DB.VOLUME; while(VOLUME != NULL) printf("%d %d\n",VOLUME->type,VOLUME->curved), VOLUME = VOLUME->next;

}
