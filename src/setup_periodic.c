// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_periodic.h"

/*
 *	Purpose:
 *		Modify connectivity if using a periodic mesh.
 *
 *	Comments:
 *		The function achieves the result using the following steps:
 *			1) Find the number of periodic facets in each direction and store their associated number of vertices so
 *			   that TRIs and QUADs are not matched up.
 *			2) The potential connections for each node are then computed and stored.
 *			3) Looping over the master facets in each direction, master vertices are compared to slave vertices on each
 *			   possible slave periodic facet using the potential connections array. If a match is found, the slave facet
 *			   is no longer considered for subsequent master facets. Once the loop has terminated, a check is performed
 *			   to ensure that the number of matches found is equal to the number of master periodic facets in the mesh
 *			   file.
 *			4) Once all periodic connections have been established, VToV, VToF, VToGF, GFToVe, and GFC are modified to
 *			   reflect the changes.
 *
 *	Notation:
 *		NPF[dim]  : (N)umber of (P)eriodic (F)acets of each (dim)ension
 *		PFNVe[]   : (P)eriodic (F)acet (N)umber of (Ve)rtices
 *		PFToGF[]  : (P)eriodic (F)acet To (G)lobal (F)acet correspondence
 *		PFVFInd[] : (P)eriodic (F)acet (V)olume and (F)acet (Ind)ices in VTo() Arrays
 *		PFToVe[]  : (P)eriodic (F)acet To (Ve)rtex correspondence
 *		ve[]      : (ve)rtices of the periodic facet under consideration
 *		v[]       : (v)olume
 *		f[]       : (f)acet
 *		            [] : (M)aster, (S)lave
 *
 *	References:
 *
 */

void setup_periodic()
{
	// Initialize DB Parameters
	unsigned int d        = DB.d,
	             NV       = DB.NV,
	             NGF      = DB.NGF,
	             NVe      = DB.NVe,
	             NPVe     = DB.NPVe,
	             NfMax    = DB.NfMax,
	             NfveMax  = DB.NfveMax,
	             *NE      = DB.NE,
	             *PVe     = DB.PVe,
	             *ETags   = DB.ETags,
	             *EType   = DB.EType,
	             *VToGF   = DB.VToGF,
	             *VToBC   = DB.VToBC,
	             *GFToVe  = DB.GFToVe,
	             *VToV    = DB.VToV,
	             *VToF    = DB.VToF,
	             *GFC     = DB.GFC;

	unsigned int PrintTesting = 0;

	// Standard datatypes
	unsigned int i, j, dim, BlockStart, count, gf, iMax,
	             IndPF, IndPFM, IndPFS, IndU, IndPVe, IndveM, IndveS, v, f, NGFC, NGFCUnique,
	             Fs, Vs, E, NPFTotal, NPVeUnique, PFM, PFS, NpveM, tmp, match, vM, vS, fM, fS,
	             NPF[d], NPFSum[d],
	             *PFTypeMOver, *PFToVeM, *PFTypeSOver, *PFToVeS, *PFNveM, *PFNveS, *PFToGFM, *PFToGFS,
	             *PFVFIndM, *PFVFIndS, *PVeUnique, *PConn, *NConn, *PFSFound, *veM, *veS, *matchM, *matchS,
	             *GFToRemove, *GFCOver, *GFToVeUnder;

	struct S_ELEMENT *ELEMENT;


	Fs = 0; for (i = 0; i < d-1; i++) Fs += NE[i];
	Vs = 0; for (i = 0; i < d; i++)   Vs += NE[i];

	for (i = 0; i < d; i++)
		NPF[i] = 0;

	PFTypeMOver = malloc(NE[d-1] * sizeof *PFTypeMOver); // free
	PFTypeSOver = malloc(NE[d-1] * sizeof *PFTypeSOver); // free

	for (dim = 0, IndPFM = 0, IndPFS = 0; dim < d; dim++) {
		for (E = Fs; E < Vs; E++) {
			if (ETags[E*2+0] % BC_STEP_SC == (PERIODIC_XL+2*dim)) {
				NPF[dim] += 1;
				PFTypeMOver[IndPFM] = EType[E];
				IndPFM++;
			} else if (ETags[E*2+0] % BC_STEP_SC == (PERIODIC_XR+2*dim)) {
				PFTypeSOver[IndPFS] = EType[E];
				IndPFS++;
			}
		}
	}
	if (IndPFM != IndPFS)
		printf("Error: Master and Slave periodic facet mismatch\n"), exit(1);

	NPFTotal = IndPFM;

	PFNveM = malloc(NPFTotal * sizeof *PFNveM); // free
	PFNveS = malloc(NPFTotal * sizeof *PFNveS); // free

	for (i = 0; i < NPFTotal; i++) {
		ELEMENT = get_ELEMENT_type(PFTypeMOver[i]);
		PFNveM[i] = ELEMENT->Nve;

		ELEMENT = get_ELEMENT_type(PFTypeSOver[i]);
		PFNveS[i] = ELEMENT->Nve;
	}
	free(PFTypeMOver);
	free(PFTypeSOver);

	for (i = 0; i < d; i++) {
		NPFSum[i] = 0;
		for (j = 0; j < i+1; j++)
			NPFSum[i] += NPF[j];
	}

//array_print_ui(1,d,NPF,'R');
//array_print_ui(1,d,NPFSum,'R');
//array_print_ui(DB.NETotal,2,ETags,'R');

	GFToRemove = malloc(NPFTotal * sizeof *GFToRemove); // free
	PFToGFM    = malloc(NPFTotal * sizeof *PFToGFM); // free
	PFToGFS    = malloc(NPFTotal * sizeof *PFToGFS); // free

	PFVFIndM = malloc(NPFTotal*2 *sizeof *PFVFIndM); // free
	PFVFIndS = malloc(NPFTotal*2 *sizeof *PFVFIndS); // free

	for (dim = 0, IndPFM = 0, IndPFS = 0; dim < d; dim++) {
		for (i = 0; i < NV; i++) {
		for (j = 0; j < NfMax; j++) {
			if (VToBC[i*NfMax+j] % BC_STEP_SC == (PERIODIC_XL+2*dim)) {
				PFVFIndM[IndPFM*2  ] = i;
				PFVFIndM[IndPFM*2+1] = j;
				PFToGFM[IndPFM]  = VToGF[i*NfMax+j];
				IndPFM++;
			} else if (VToBC[i*NfMax+j] % BC_STEP_SC == (PERIODIC_XR+2*dim)) {
				PFVFIndS[IndPFS*2  ] = i;
				PFVFIndS[IndPFS*2+1] = j;
				PFToGFS[IndPFS]  = VToGF[i*NfMax+j];
				IndPFS++;
			}
		}}
	}

//array_print_ui(1,NPFTotal,PFToGFM,'R');
//array_print_ui(1,NPFTotal,PFToGFS,'R');

	PFToVeM = malloc(NPFSum[d-1]*NfveMax * sizeof *PFToVeM); // free
	PFToVeS = malloc(NPFSum[d-1]*NfveMax * sizeof *PFToVeS); // free
	for (i = 0; i < NPFSum[d-1]*NfveMax; i++) {
		PFToVeM[i] = NVe;
		PFToVeS[i] = NVe;
	}

	for (dim = 0, IndPF = 0, BlockStart = 0; dim < d; dim++) {
		if (dim != 0)
			BlockStart += NPF[dim-1]*NfveMax;
		for (i = 0; i < NPF[dim]; i++) {
			for (j = 0; j < PFNveM[IndPF]; j++)
				PFToVeM[BlockStart+i*NfveMax+j] = GFToVe[PFToGFM[IndPF]*NfveMax+j];
			PetscSortInt(PFNveM[IndPF],(int *)&PFToVeM[BlockStart+i*NfveMax]);

			for (j = 0; j < PFNveS[IndPF]; j++)
				PFToVeS[BlockStart+i*NfveMax+j] = GFToVe[PFToGFS[IndPF]*NfveMax+j];
			PetscSortInt(PFNveS[IndPF],(int *)&PFToVeS[BlockStart+i*NfveMax]);

			IndPF++;
		}
	}
	free(PFToGFM);
	free(PFToGFS);

//array_print_ui(NPFSum[d-1],NfveMax,PFToVeM,'R');
//array_print_ui(NPFSum[d-1],NfveMax,PFToVeS,'R');
//array_print_ui(DB.NPVe,2,DB.PVe,'R');

	// Make array of possible connections for each node
	for (i = NPVeUnique = 1, IndU = 0; i < NPVe; i++) {
		if (PVe[i*2+0] != PVe[IndU*2+0]) {
			NPVeUnique++;
			IndU = i;
		}
	}

	PVeUnique = malloc(NPVeUnique   * sizeof *PVeUnique); // free
	PConn     = malloc(NPVeUnique*3 * sizeof *PConn); // free
	NConn     = malloc(NPVeUnique*3 * sizeof *NConn); // free
	for (i = 0, iMax = NPVeUnique*3; i < iMax; i++) PConn[i] = NVe;

	for (i = IndPVe = 0; i < NPVeUnique; i++) {
		PVeUnique[i] = PVe[IndPVe*2+0];
		j = 0;
		do {
			PConn[i*3+j] = PVe[(IndPVe+j)*2+1];
			j++;
		} while((IndPVe+j) < NPVe && PVe[IndPVe*2+0] == PVe[(IndPVe+j)*2+0]);
		NConn[i] = j;
		IndPVe += j;
	}

//array_print_ui(1,NPVeUnique,NConn,'R');
//array_print_ui(NPVeUnique,3,PConn,'R');
//array_print_ui(1,NPVeUnique,PVeUnique,'R');

	// Find Corresponding Periodic Facets
	for (dim = 0, BlockStart = 0, IndPFM = 0; dim < d; dim++) {
		if (dim != 0)
			BlockStart += NPF[dim-1]*NfveMax;

		PFSFound = calloc(NPF[dim] , sizeof *PFSFound); // free
		for (PFM = 0; PFM < NPF[dim]; PFM++) {
			// For each periodic master facet, find the corresponding slave facet
			NpveM = PFNveM[IndPFM];
			veM = malloc(NpveM * sizeof *veM); // free
			for (IndveM = 0; IndveM < NpveM; IndveM++)
				veM[IndveM] = PFToVeM[IndPFM*NfveMax+IndveM];

			veS = malloc(NpveM * sizeof *veS); // free
			for (i = IndPFS = 0; i < dim; i++)
				IndPFS += NPF[i];

			for (PFS = 0; PFS < NPF[dim]; PFS++) {
				if (PFNveS[IndPFS] != NpveM)
					continue; // ensure that TRIs and QUADs are not paired.

				for (i = 0; i < NpveM; i++)
					veS[i] = PFToVeS[IndPFS*NfveMax+i];

				matchM = malloc(NpveM * sizeof *matchM); // free
				matchS = malloc(NpveM * sizeof *matchS); // free
				for (i = 0; i < NpveM; i++) {
					matchM[i] = NVe;
					matchS[i] = NVe;
				}

				// ToBeDeleted: Need to make sure that this never links the wrong facets (THINK)
				for (IndveM = 0; IndveM < NpveM; IndveM++) {
				for (IndveS = 0; IndveS < NpveM; IndveS++) {
					if (PFSFound[PFS] == 0 && matchM[IndveM] == NVe && matchS[IndveS] == NVe) {
						array_find_indexo_ui(NPVeUnique,PVeUnique,veS[IndveS],&IndPVe,&tmp);
						for (i = 0; i < NConn[IndPVe]; i++) {
							if (PConn[IndPVe*3+i] == veM[IndveM]) {
								matchM[IndveM] = veM[IndveM];
								matchS[IndveS] = veS[IndveS];
								break;
							}
						}
					}
				}}

				for (i = 0, match = 1; i < NpveM; i++) {
					if (matchS[i] == NVe) {
						match = 0;
						break;
					}
				}
				free(matchM);
				free(matchS);

				// Found the matching periodic facet
				if (match) {
					PFSFound[PFS] = 1;
					// Replace slave facet information in VToV and VToGF with master facet information
					vM = PFVFIndM[IndPFM*2];
					fM = PFVFIndM[IndPFM*2+1];
					vS = PFVFIndS[IndPFS*2];
					fS = PFVFIndS[IndPFS*2+1];

//printf("%d %d %d %d \n",vM,fM,vS,fS);

					VToV[vM*NfMax+fM] = vS;
					VToF[vM*NfMax+fM] = fS;
					VToV[vS*NfMax+fS] = vM;
					VToF[vS*NfMax+fS] = fM;

					GFToRemove[IndPFM] = VToGF[vS*NfMax+fS];
					VToGF[vS*NfMax+fS] = VToGF[vM*NfMax+fM];

					IndPFM++;
					break;
				}
				IndPFS++;

			}
			free(veM);
			free(veS);
		}
		free(PFSFound);
	}
	free(PFNveM);
	free(PFNveS);
	free(PFVFIndM);
	free(PFVFIndS);
	free(PFToVeM);
	free(PFToVeS);
	free(PVeUnique);
	free(PConn);
	free(NConn);

	if (IndPFM != NPFSum[d-1])
		printf("Error: Did not find the correct number of periodic connections.\n"), exit(1);

	PetscSortInt(NPFSum[d-1],(int *)GFToRemove);

	// Decrement VToGF
	for (v = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		count = 0;
		while(count < NPFSum[d-1] && VToGF[v*NfMax+f] > GFToRemove[count])
			count++;
		VToGF[v*NfMax+f] -= count;
	}}

	// Fix and reallocate memory for other affected arrays
	GFCOver = malloc(NV*NfMax * sizeof *GFCOver); // free
	for (v = 0, NGFC = 0; v < NV; v++) {
	for (f = 0; f < NfMax; f++) {
		if (VToBC[v*NfMax+f] >= 2000) {
			GFCOver[NGFC] = VToGF[v*NfMax+f];
			NGFC++;
		}
	}}
	PetscSortInt(NGFC,(int *)GFCOver);

	NGFCUnique = 1;
	for (i = 1; i < NGFC; i++) {
		if (GFCOver[i] != GFCOver[NGFCUnique-1]) {
			GFCOver[NGFCUnique] = GFCOver[i];
			NGFCUnique++;
		}
	}

	free(GFC);
	GFC = malloc(NGFCUnique * sizeof *GFC); // keep
	for (i = 0; i < NGFCUnique; i++)
		GFC[i] = GFCOver[i];
	free(GFCOver);

	NGF -= NPFSum[d-1];
	GFToVeUnder = malloc(NGF*NfveMax * sizeof *GFToVeUnder); // keep
	for (gf = 0, count = 0; gf < NGF; gf++) {
		while (gf+count == GFToRemove[count])
			count++;
		for (j = 0; j < NfveMax; j++)
			GFToVeUnder[gf*NfveMax+j] = GFToVe[(gf+count)*NfveMax+j];
	}
	free(GFToVe);

	DB.NGF    = NGF;
	DB.NGFC   = NGFCUnique;
	DB.GFC    = GFC;
	DB.GFToVe = GFToVeUnder;

	// Testing
	if (PrintTesting && DB.Testing && !DB.MPIrank) {
		printf("VToV:\n");       array_print_ui(DB.NV,DB.NfMax,DB.VToV,'R');
		printf("VToF:\n");       array_print_ui(DB.NV,DB.NfMax,DB.VToF,'R');
		printf("VToGF:\n");      array_print_ui(DB.NV,DB.NfMax,DB.VToGF,'R');
		//printf("GFToRemove:\n"); array_print_ui(1,NPFSum[d-1],GFToRemove,'R');
		printf("GFToVe:\n");     array_print_ui(DB.NGF,DB.NfveMax,DB.GFToVe,'R');
		printf("GFC:\n");            array_print_ui(1,DB.NGFC,GFC,'R');
	}

	free(GFToRemove);
	//free(DB.ETags); // No longer needed, but this should not be freed only if periodic (ToBeModified)
}
