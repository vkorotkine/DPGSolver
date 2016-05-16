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
 *		Will need two different setup_structures functions: One using the initial global arrays and one updating
 *		elements individually after hp refinement. (ToBeDeleted)
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
				 NfMax   = DB.NfMax,
	             NV      = DB.NV,
	             NGF     = DB.NGF,
	             NGFC    = DB.NGFC,
	             NVC     = DB.NVC,
	             *NE     = DB.NE,
	             *EToVe  = DB.EToVe,
	             *EType  = DB.EType,
	             *EToPrt = DB.EToPrt,
//	             *VToV   = DB.VToV,
	             *VToGF  = DB.VToGF,
	             *VC     = DB.VC,
	             *GFC    = DB.GFC;
	int          MPIrank = DB.MPIrank;
	double       *VeXYZ  = DB.VeXYZ;

	int  PrintTesting = 0;

	// Standard datatypes
	unsigned int i, iMax, f, v, dim, ve, gf, curved,
	             IndE, IndVC, IndVgrp, IndGFC,
	             Vs, vlocal, NVlocal, NECgrp, NVgrp,
				 indexg, NvnGs,
	             uMPIrank,
	             *GFToV, *GF_Nv;
	double       *XYZc;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **Vgrp, **Vgrp_tmp;
	struct S_FACET   **FACET, **FoundFACET;

	// silence
	NECgrp = 0;

	uMPIrank = MPIrank;

	if      (d == 1) NECgrp = 1;
	else if (d == 2) NECgrp = 2;
	else if (d == 3) NECgrp = 4;

	FACET      = calloc(2   , sizeof *FACET); // free
	FoundFACET = calloc(NGF , sizeof *FoundFACET); // free

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
	DB.FACET  = NULL;

	// Note: only initialize volumes on the current processor.
	for (v = 0, NVlocal = 0; v < NV; v++) {
		if (EToPrt[v] == uMPIrank)
			NVlocal++;
	}

	for (v = 0, IndE = Vs, IndGFC = IndVC = 0, vlocal = 0; v < NV; v++) {
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

			// FACETs adjacent to VOLUMEs on the current processor.
			// Note that GFC and VToGF cycle through the global FACET indices in order
			for (f = 0; f < NfMax; f++) {
				gf = VToGF[v*NfMax+f];
				if (FoundFACET[gf] == NULL) {

					if (DB.FACET != NULL) {
						FACET[0]->next = New_FACET();
						FACET[0]       = FACET[0]->next;
					} else {
						DB.FACET = New_FACET();
						FACET[0] = DB.FACET;
					}

//					FACET[0]->indexl = gflocal;
					FACET[0]->indexg = gf;
					FACET[0]->P      = VOLUME->P;

					FACET[0]->VIn   = VOLUME;
					FACET[0]->VfIn  = 9*f;

					// Overwritten if a second VOLUME is found adjacent to this FACET
					FACET[0]->VOut  = VOLUME;
					FACET[0]->VfOut = 9*f;

					if (!VOLUME->curved) {
						FACET[0]->typeInt = 's';
					} else {
						FACET[0]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACET[0]->curved = 1;
							IndGFC++;
						}
					}

					FoundFACET[gf] = FACET[0];
				} else {
					FACET[1] = FoundFACET[gf];

					FACET[1]->P = max(FACET[1]->P,VOLUME->P);
					FACET[1]->VOut  = VOLUME;
					FACET[1]->VfOut = 9*f;
					if (VOLUME->curved) {
						FACET[1]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACET[1]->curved = 1;
							IndGFC++;
						}
					}
				}
//				// Indexing from connectivity discussion above noting that the mesh is conforming at the start (ToBeDeleted)
//				VOLUME->GF[f*9] = FoundFACET[gf];
			}


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
			for (f = 0; f < NfMax; f++) {
				gf = VToGF[v*NfMax+f];
				if (IndGFC < NGFC && gf == GFC[IndGFC])
					IndGFC++;
			}
		}

		IndE++;
	}
	free(Vgrp_tmp);
	free(FACET);
	free(FoundFACET);

	if (!AC && IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

/*
	// Determine GFToV array
	GFToV = malloc(NGF*2 * sizeof *GFToV); // tbd
	GF_Nv = calloc(NGF   , sizeof *GF_Nv); // free

	for (gf = 0; gf < NGF; gf++)
		GFToV[gf*2+1] = NV;

	for (v = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		gf = VToGF[v*NfMax+f];

		GFToV[gf*2+GF_Nv[gf]] = v;
		GF_Nv[gf]++;
	}}
	free(GF_Nv);

array_print_ui(NGF,2,GFToV,'R');
exit(1);
*/

/*
for (FACET[0] = DB.FACET; FACET[0] != NULL; FACET[0] = FACET[0]->next) {
	printf("%d %d %c %d %d %d %d\n",
	       FACET[0]->indexg,FACET[0]->curved,FACET[0]->typeInt,FACET[0]->VIn->indexg,FACET[0]->VOut->indexg,FACET[0]->VfIn,FACET[0]->VfOut);
}
exit(1);
*/


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
