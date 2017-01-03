// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_structures.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
 
#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "array_norm.h"
#include "memory_constructors.h"
#include "element_functions.h"
#include "matrix_functions.h"

/*
 *	Purpose:
 *		Set up VOLUME and FACE structures.
 *
 *	Comments:
 *		Will need two different setup_structures functions: One using the initial global arrays and one updating
 *		elements individually after hp refinement. (ToBeDeleted)
 *
 *	Notation:
 *		XYZ_(1)(2) : Physical node locations (XYZ) of (1) nodes of (2) type
 *		             (1) : (v)olume, (f)ace
 *		             (2) : (V)ertex
 *		IndOrd(In/Out)(Out/In) : (Ind)ex of (Ord)ering relating projection from (In)ner VOLUME to the FACE to the
 *		                         projection from the (Out)er VOLUME to the FACE
 *
 *	References:
*/

static void compute_distance_matrix(const unsigned int Nn, const unsigned int BC, const unsigned int d, double *XYZIn,
                                    double *XYZOut, double *DXYZ)
{
	// Standard datatypes
	unsigned int i, j, k, kMax,
	             tmp_ui, IndDXYZ, IndComp[2];
	double       tmp_d;

	// silence
	IndComp[0] = IndComp[1] = 0;

	if (d == 1)
		return;

	for (i = 0; i < Nn; i++) {
	for (j = 0; j < Nn; j++) {
		tmp_ui = BC % BC_STEP_SC;
		if (tmp_ui > BC_PERIODIC_MIN) { // special case for periodic
			if (d == 3) {
				if      (tmp_ui == PERIODIC_XL || tmp_ui == PERIODIC_XR) IndComp[0] = 1, IndComp[1] = 2;
				else if (tmp_ui == PERIODIC_YL || tmp_ui == PERIODIC_YR) IndComp[0] = 0, IndComp[1] = 2;
				else if (tmp_ui == PERIODIC_ZL || tmp_ui == PERIODIC_ZR) IndComp[0] = 0, IndComp[1] = 1;
			} else if (d == 2) {
				if      (tmp_ui == PERIODIC_XL || tmp_ui == PERIODIC_XR) IndComp[0] = 1;
				else if (tmp_ui == PERIODIC_YL || tmp_ui == PERIODIC_YR) IndComp[0] = 0;
			}
			IndDXYZ = i*Nn+j;
			for (k = 0, kMax = d-1; k < kMax; k++) {
				tmp_d = fabs(XYZIn[IndComp[k]*Nn+i]-XYZOut[IndComp[k]*Nn+j]);
				if (tmp_d > DXYZ[IndDXYZ])
					DXYZ[IndDXYZ] = tmp_d;
			}
		} else {
			IndDXYZ = i*Nn+j;
			for (k = 0; k < d; k++)
				DXYZ[IndDXYZ] += fabs(XYZIn[k*Nn+i]-XYZOut[k*Nn+j]);
		}
	}}
}

static void get_ordering_index(const unsigned int Nn, const unsigned int d, double *DXYZ, unsigned int *IndOrdInOut,
                               unsigned int *IndOrdOutIn)
{
	/*
	 *	Purpose:
	 *		Return the ordering index corresponding to the match between two FACEs based on DXYZ.
	 *
	 *	Comments:
	 *		In 1D, IndOrdInOut == IndOrdOutIn == 0.
	 *		In 2D, IndOrdInOut == IndOrdOutIn.
	 *
	 *	Notation:
	 *		IndZerosP : (Ind)ices of (Zeros) which are (P)ossible.
	 */

	if (d == 1) {
		*IndOrdInOut = 0;
		*IndOrdOutIn = 0;
		return;
	} else {
		unsigned int i, j, countZeros,
					 IndZerosInOut[Nn], IndZerosOutIn[Nn];

		// Find indices of zeros in DXYZ
		countZeros = 0;
		for (i = 0; i < Nn; i++) {
		for (j = 0; j < Nn; j++) {
			if (fabs(DXYZ[i*Nn+j]) < EPS) {
				IndZerosOutIn[i] = j;
				countZeros++;
				break;
			}
		}}
		if (countZeros != Nn)
			printf("Error: Did not find a sufficient number of zeros in DXYZ (get_ordering_index).\n"), exit(1);

		if (d == 3) {
			// Transpose DXYZ and find Out->In Ordering as well
			mkl_dimatcopy('R','T',Nn,Nn,1.0,DXYZ,Nn,Nn);

			// Find indices of zeros in DXYZ'
			for (i = 0; i < Nn; i++) {
			for (j = 0; j < Nn; j++) {
				if (fabs(DXYZ[i*Nn+j]) < EPS) {
					IndZerosInOut[i] = j;
					break;
				}
			}}

			if (Nn == 4) { // QUAD FACE
				unsigned int IndZerosP[32] = { 0, 1, 2, 3,
				                               1, 0, 3, 2,
				                               2, 3, 0, 1,
				                               3, 2, 1, 0,
				                               0, 2, 1, 3,
				                               2, 0, 3, 1,
				                               1, 3, 0, 2,
				                               3, 1, 2, 0};

				for (i = 0; i < 8; i++) {
					if (array_norm_diff_ui(Nn,IndZerosInOut,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdInOut = i;
					if (array_norm_diff_ui(Nn,IndZerosOutIn,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdOutIn = i;
				}
			} else if (Nn == 3) { // TRI FACE
				unsigned int IndZerosP[18] = { 0, 1, 2,
				                               1, 2, 0,
				                               2, 0, 1,
				                               0, 2, 1,
				                               2, 1, 0,
				                               1, 0, 2};

				for (i = 0; i < 6; i++) {
					if (array_norm_diff_ui(Nn,IndZerosInOut,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdInOut = i;
					if (array_norm_diff_ui(Nn,IndZerosOutIn,&IndZerosP[i*Nn],"Inf") < EPS)
						*IndOrdOutIn = i;
				}
			}
		} else if (d == 2) {
			unsigned int IndZerosP[4] = { 0, 1,
			                              1, 0};
			for (i = 0; i < 2; i++) {
				if (array_norm_diff_ui(Nn,IndZerosOutIn,&IndZerosP[i*Nn],"Inf") < EPS) {
					*IndOrdInOut = i;
					*IndOrdOutIn = i;
					return;
				}
			}
		}
	}
}

void setup_structures(void)
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             NfrefMax = DB.NfrefMax,
	             AC       = DB.AC,
	             PGlobal  = DB.PGlobal,
	             NP       = DB.NP,
	             NfMax    = DB.NfMax,
	             NV       = DB.NV,
	             NGF      = DB.NGF,
	             NGFC     = DB.NGFC,
	             NVC      = DB.NVC,
	             *NE      = DB.NE,
	             *EToVe   = DB.EToVe,
	             *EType   = DB.EType,
	             *EToPrt  = DB.EToPrt,
	             *VToV    = DB.VToV,
	             *VToF    = DB.VToF,
	             *VToGF   = DB.VToGF,
	             *VToBC   = DB.VToBC,
	             *VC      = DB.VC,
	             *GFC     = DB.GFC;
	int          MPIrank  = DB.MPIrank;
	double       *VeXYZ   = DB.VeXYZ;

	// Standard datatypes
	unsigned int i, f, v, dim, ve, gf, fh, fhMax,
	             IndE, IndVC, IndVgrp, IndGFC, IndVIn, Indf, IndOrdInOut, IndOrdOutIn,
	             Vs, vlocal, NTVgrp, NVlocal, NECgrp, *NVgrp, Vf, Nf,
	             VfIn, fIn, VfOut, fOut,
	             Nve,*Nfve, Nfn,
				 indexg, NvnGs,
	             uMPIrank;
	double       *XYZ_vV, **VeF, *XYZIn_fC, *XYZOut_fC, *DXYZ;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME, **Vgrp, **Vgrp_tmp, *VIn, *VOut;
	struct S_FACE   **FACE, **FoundFACE;

	// silence
	Nf = NECgrp = IndOrdInOut = IndOrdOutIn = 0;

	uMPIrank = MPIrank;

	if      (d == 1) NECgrp = 1;
	else if (d == 2) NECgrp = 2;
	else if (d == 3) NECgrp = 4;

	FACE      = calloc(2   , sizeof *FACE); // free
	FoundFACE = calloc(NGF , sizeof *FoundFACE); // free

	NTVgrp = NECgrp*NP*2;
	NVgrp    = calloc(NTVgrp , sizeof *NVgrp);    // keep
	Vgrp     = calloc(NTVgrp , sizeof *Vgrp);     // keep
	Vgrp_tmp = calloc(NTVgrp , sizeof *Vgrp_tmp); // free

	Vs = 0; for (i = 0; i < d; i++) Vs += NE[i];

	VOLUME = New_VOLUME();
	DB.VOLUME = VOLUME;
	DB.FACE  = NULL;

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
			VOLUME->P      = PGlobal;
			VOLUME->PNew   = PGlobal;
			VOLUME->level  = 0;
			VOLUME->type   = EType[IndE];
			VOLUME->Eclass = get_Eclass(VOLUME->type);
			VOLUME->update = 1;

			if (AC || (IndVC < NVC && v == VC[IndVC])) {
				VOLUME->curved = 1;
				IndVC++;
			} else {
				VOLUME->curved = 0;
			}

			// FACEs adjacent to VOLUMEs on the current processor.
			// Note that GFC and VToGF cycle through the global FACE indices in order
			ELEMENT = get_ELEMENT_type(VOLUME->type);

			Nf = ELEMENT->Nf;

			for (f = 0; f < Nf; f++) {
				VOLUME->NsubF[f] = 1;
				gf = VToGF[v*NfMax+f];
				if (FoundFACE[gf] == NULL) {
					if (DB.FACE) {
						FACE[0]->next = New_FACE();
						FACE[0]       = FACE[0]->next;
					} else {
						DB.FACE = New_FACE();
						FACE[0] = DB.FACE;
					}

//					FACE[0]->indexl = gflocal;
					FACE[0]->indexg = gf;
					FACE[0]->P      = VOLUME->P;
					FACE[0]->level  = 0;

					FACE[0]->VIn   = VOLUME; IndVIn = VOLUME->indexg;
					FACE[0]->VfIn  = NfrefMax*f;

					FACE[0]->BC    = VToBC[IndVIn*NfMax+f];

					// Overwritten if a second VOLUME is found adjacent to this FACE
					FACE[0]->VOut  = VOLUME;
					FACE[0]->VfOut = NfrefMax*f;

					if (!VOLUME->curved) {
						FACE[0]->typeInt = 's';
					} else {
						FACE[0]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACE[0]->curved = 1;
							IndGFC++;
						}
					}

					FACE[0]->VOut->FACE[f*NSUBFMAX] = FACE[0];

					FoundFACE[gf] = FACE[0];
				} else {
					FACE[1] = FoundFACE[gf];

					FACE[1]->P = max(FACE[1]->P,VOLUME->P);
					FACE[1]->VOut  = VOLUME;
					FACE[1]->VfOut = NfrefMax*f;
					FACE[1]->VOut->FACE[f*NSUBFMAX] = FACE[1];
					if (VOLUME->curved) {
						FACE[1]->typeInt = 'c';
						if (AC || (IndGFC < NGFC && gf == GFC[IndGFC])) {
							FACE[1]->curved = 1;
							IndGFC++;
						}
					}
				}
//				// Indexing from connectivity discussion above noting that the mesh is conforming at the start (ToBeDeleted)
//				VOLUME->GF[f*NfrefMax] = FoundFACE[gf];
			}


			// Geometry
			if (VOLUME->Eclass == C_TP)
				NvnGs = pow(ELEMENT->ELEMENTclass[0]->NvnGs[1],d);
			else if (VOLUME->Eclass == C_WEDGE)
				NvnGs = (ELEMENT->ELEMENTclass[0]->NvnGs[1])*(ELEMENT->ELEMENTclass[1]->NvnGs[1]);
			else if (VOLUME->Eclass == C_SI || VOLUME->Eclass == C_PYR)
				NvnGs = ELEMENT->NvnGs[1];
			else
				printf("Error: Unsupported element type setup_struct (NvnGs).\n"), exit(1);

			XYZ_vV = malloc(NvnGs*d * sizeof *XYZ_vV); // keep
			VOLUME->XYZ_vV = XYZ_vV;

			// XYZ_vV may be interpreted as [X Y Z] where each of X, Y, Z are column vectors (ToBeDeleted)
			// Add this comment to notation section.
			indexg = VOLUME->indexg;
			for (ve = 0; ve < NvnGs; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZ_vV[dim*NvnGs+ve] = VeXYZ[EToVe[(Vs+indexg)*NVEMAX+ve]*d+dim];
			}}
//printf("%d\n",VOLUME->indexg);
//array_print_d(NvnGs,d,XYZ_vV,'C');

			// MPI
			IndVgrp = ((VOLUME->Eclass)*NP*2)+(VOLUME->P*2)+(VOLUME->curved);
			if (!NVgrp[IndVgrp])
				Vgrp[IndVgrp] = VOLUME;
			else
				Vgrp_tmp[IndVgrp]->grpnext = VOLUME;

			NVgrp[IndVgrp]++;
			Vgrp_tmp[IndVgrp] = VOLUME;

			if (vlocal != NVlocal-1) {
				VOLUME->next = New_VOLUME();
				VOLUME = VOLUME->next;
			}

			vlocal++;
		} else {
			// Ensure that appropriate global indices are incremented if necessary
			if (NVC && v == VC[IndVC])
				IndVC++;
			for (f = 0; f < Nf; f++) {
				gf = VToGF[v*NfMax+f];
				if (IndGFC < NGFC && gf == GFC[IndGFC])
					IndGFC++;
			}
		}

		IndE++;
	}
	free(Vgrp_tmp);
	free(FoundFACE);

	if (!AC && IndVC > NVC)
		printf("Error: Found too many curved VOLUMEs.\n"), exit(1);

	// Flag boundary FACEs
	for (FACE[0] = DB.FACE; FACE[0]; FACE[0] = FACE[0]->next) {
		VIn   = FACE[0]->VIn;
		VfIn  = FACE[0]->VfIn;
		fIn   = VfIn/NFREFMAX;

		VOut  = FACE[0]->VOut;
		VfOut = FACE[0]->VfOut;
		fOut  = VfOut/NFREFMAX;

		FACE[0]->Boundary = !((VIn->indexg != VOut->indexg) || (VIn->indexg == VOut->indexg && fIn != fOut));
	}

	// Initialize VOLUME connectivity
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		// Initialize to NV (impossible value)
		for (fh = 0, fhMax = NFMAX*NFREFMAX; fh < fhMax; fh++)
			VOLUME->neigh[fh] = NV;

		indexg = VOLUME->indexg;
		ELEMENT = get_ELEMENT_type(VOLUME->type);
		Nf = ELEMENT->Nf;
		for (f = 0; f < Nf; f++) {
			VOLUME->neigh[f*NFREFMAX]   = VToV[NfMax*indexg+f];
			VOLUME->neigh_f[f*NFREFMAX] = VToF[NfMax*indexg+f];
		}
	}

	for (FACE[0] = DB.FACE; FACE[0]; FACE[0] = FACE[0]->next) {
// May potentially have a problem for PYR-HEX interface due to PYR nodes being ordered for symmetry while QUAD nodes are
// ordered for TP extension. (ToBeDeleted)
// Probably not as the FACE nodes would not be related to the PYR node ordering.

		// Obtain FACE type
		VOLUME  = FACE[0]->VIn;
		Vf      = FACE[0]->VfIn;
		Indf    = Vf / NfrefMax; // face index (ToBeDeleted: Move to notation)

		if (d == 1) {
			FACE[0]->type = POINT;
		} else if (d == 2) {
			FACE[0]->type = LINE;
		} else if (d == 3) {
			VIn  = FACE[0]->VIn;

			if (VIn->type == TET || (VIn->type == WEDGE && Indf > 2) || (VIn->type == PYR && Indf < 4))
				FACE[0]->type = TRI;
			else
				FACE[0]->type = QUAD;
		}

		// Obtain XYZIn_fC/XYZOut_fC
		VOLUME  = FACE[0]->VIn;
		Vf      = FACE[0]->VfIn;
		Indf    = Vf / NfrefMax; // face index (ToBeDeleted: Move to notation)

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nve   = ELEMENT->Nve;
		Nfve  = ELEMENT->Nfve;
		VeF   = ELEMENT->VeF;

		NvnGs = ELEMENT->NvnGs[1];
		XYZ_vV = VOLUME->XYZ_vV;

		XYZIn_fC = malloc(Nfve[Indf]*d * sizeof *XYZIn_fC); // free
		mm_CTN_d(Nfve[Indf],d,Nve,VeF[Vf],XYZ_vV,XYZIn_fC);

		VOLUME  = FACE[0]->VOut;
		Vf      = FACE[0]->VfOut;
		Indf    = Vf / NfrefMax; // face index (ToBeDeleted: Move to notation)

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nve   = ELEMENT->Nve;
		Nfve  = ELEMENT->Nfve;
		VeF   = ELEMENT->VeF;

		NvnGs = ELEMENT->NvnGs[1];
		XYZ_vV = VOLUME->XYZ_vV;

		XYZOut_fC = malloc(Nfve[Indf]*d * sizeof *XYZOut_fC); // free
		mm_CTN_d(Nfve[Indf],d,Nve,VeF[Vf],XYZ_vV,XYZOut_fC);

		// Compute distance matrix
		Nfn  = Nfve[Indf];
		DXYZ = calloc(Nfn*Nfn , sizeof *DXYZ); // free
		compute_distance_matrix(Nfn,FACE[0]->BC,d,XYZIn_fC,XYZOut_fC,DXYZ);

		// Obtain the index of corresponding ordering between FACEs
		get_ordering_index(Nfn,d,DXYZ,&IndOrdInOut,&IndOrdOutIn);

		FACE[0]->IndOrdInOut = IndOrdInOut;
		FACE[0]->IndOrdOutIn = IndOrdOutIn;

/*
printf("\n\n%d\n",FACE[0]->indexg);
array_print_d(Nfve[Indf],d,XYZIn_fC,'C');
array_print_d(Nfve[Indf],d,XYZOut_fC,'C');

array_print_d(Nfn,Nfn,DXYZ,'R');

printf("InOut: %d %d\n",IndOrdInOut,IndOrdOutIn);
*/

		free(XYZIn_fC);
		free(XYZOut_fC);

		free(DXYZ);
	}
	free(FACE);

/*
for (i = 0, iMax = NTVgrp; iMax--; i++) {
	for (VOLUME = Vgrp[i]; VOLUME; VOLUME = VOLUME->grpnext) {
		printf("%d %d %d %d\n",i,VOLUME->Eclass,VOLUME->P,VOLUME->curved);
	}
	printf("\t\t%p\n",Vgrp[i]);
}
*/
//exit(1);


	// Assign/Overwrite DB parameters
	DB.NVglobal = NV;
	DB.NV       = NVlocal;
	DB.NTVgrp   = NTVgrp,
	DB.NECgrp   = NECgrp;
	DB.NVgrp    = NVgrp;
	DB.Vgrp     = Vgrp;

//VOLUME = DB.VOLUME; while(VOLUME) printf("%d %d\n",VOLUME->type,VOLUME->curved), VOLUME = VOLUME->next;
}
