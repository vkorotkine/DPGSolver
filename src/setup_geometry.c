// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_geometry.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "output_to_paraview.h"
#include "setup_ToBeCurved.h"
#include "setup_Curved.h"
#include "setup_geom_factors.h"
#include "setup_normals.h"
#include "vertices_to_exact_geom.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Set up geometry related information.
 *
 *	Comments:
 *		Loops are placed outside of the functions listed below (ToBeModified) as they are reused to update individual
 *		VOLUMEs after refinement:
 *			setup_ToBeCurved
 *			setup_geom_factors
 *		This should likely be changed (i.e. put loops inside for efficiency) as a different function will be used to
 *		updated than to set up. (ToBeDeleted)
 *
 *	Notation:
 *		XYZ_S  : High-order VOLUME (XYZ) coordinates at (S)tart (i.e. before curving)
 *		VeInfo : (Info)rmation relating to (Ve)rtices by column (#):
 *		         (0): curved
 *		         (1): update
 *		         (2): curved surface index
 *
 *	References:
*/

struct S_OPERATORS {
	unsigned int NvnG, NfnI, NfnS;
	double       **I_vG_fI, **I_vG_fS;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass);

void setup_straight(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int NvnG;
	double       *XYZ, *XYZ_S;

	NvnG  = VOLUME->NvnG;
	XYZ_S = VOLUME->XYZ_S;

	XYZ = malloc(NvnG*d * sizeof *XYZ); // keep
	for (unsigned int i = 0, iMax = d*NvnG; i < iMax; i++)
		XYZ[i] = XYZ_S[i];
	VOLUME->XYZ = XYZ;
}

void setup_FACET_XYZ(struct S_FACET *FACET)
{
	// Initialize DB Parameters
	unsigned int d     = DB.d,
	             Adapt = DB.Adapt;

	// Standard datatypes
	unsigned int VfIn, fIn, Eclass, IndFType, NvnG, NfnS, NfnI;
	double       *XYZ_fS, *XYZ_fI;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VIn;

	OPS = malloc(sizeof *OPS); // free

	VIn  = FACET->VIn;
	VfIn = FACET->VfIn;
	fIn  = VfIn/NFREFMAX;

	Eclass = get_Eclass(VIn->type);
	IndFType = get_IndFType(Eclass,fIn);

	init_ops(OPS,VIn,FACET,IndFType);

	NvnG = OPS->NvnG;
	switch (Adapt) {
	default: // ADAPT_P, ADAPT_H, ADAPT_HP
printf("Error: Should not be entering default in setup_FACET_XYZ.\n"), exit(1);
		NfnS = OPS->NfnS;

		XYZ_fS = malloc(NfnS*d *sizeof *XYZ_fS); // keep
		mm_CTN_d(NfnS,d,NvnG,OPS->I_vG_fS[VfIn],VIn->XYZ,XYZ_fS);

		FACET->XYZ_fS = XYZ_fS;
		break;
case ADAPT_P: // ToBeModified
case ADAPT_H:
case ADAPT_HP:
	case ADAPT_0:
		NfnI = OPS->NfnI;

		XYZ_fI = malloc(NfnI*d *sizeof *XYZ_fI); // keep
		mm_CTN_d(NfnI,d,NvnG,OPS->I_vG_fI[VfIn],VIn->XYZ,XYZ_fI);

		if (FACET->XYZ_fI)
			free(FACET->XYZ_fI);
		FACET->XYZ_fI = XYZ_fI;
		break;
	}
	free(OPS);
}

static void mark_curved_vertices()
{
	// Initialize DB Parameters
	unsigned int d      = DB.d,
	             NVe    = DB.NVe,
	             *NE    = DB.NE,
	             *EToVe = DB.EToVe;

	// Standard datatypes
	unsigned int dim, f, ve, Indve, Vs, *Nfve, *VeFcon, *VeInfo;

	struct S_ELEMENT *ELEMENT;
	struct S_FACET   *FACET;
	struct S_VOLUME  *VOLUME;

	Vs = 0; for (dim = 0; dim < d; dim++) Vs += NE[dim];

	VeInfo = calloc(NVe*NVEINFO , sizeof *VeInfo); // keep
	DB.VeInfo = VeInfo;

	for (ve = 0; ve < NVe; ve++)
		VeInfo[2*NVe+ve] = UINT_MAX;

	for (FACET = DB.FACET; FACET; FACET = FACET->next) {
		if (!FACET->curved)
			continue;

		VOLUME = FACET->VIn;
		f      = (FACET->VfIn)/NFREFMAX;

		ELEMENT = get_ELEMENT_type(VOLUME->type);

		Nfve   = ELEMENT->Nfve;
		VeFcon = ELEMENT->VeFcon;

		for (ve = 0; ve < Nfve[f]; ve++) {
			Indve = EToVe[(Vs+(VOLUME->indexg))*NVEMAX+VeFcon[f*NFVEMAX+ve]];
			VeInfo[0*NVe+Indve] = 1;
			VeInfo[1*NVe+Indve] = 1;
		}
	}
}

static void initialize_VOLUME_VeInfo(void)
{
	// Initialize DB Parameters
	unsigned int d       = DB.d,
	             NVe     = DB.NVe,
	             *NE     = DB.NE,
	             *VeInfo = DB.VeInfo,
	             *EToVe  = DB.EToVe;

	// Standard datatypes
	unsigned int i, dim, ve, Indve, Vs, Nve, *VOL_VeInfo;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	Vs = 0; for (dim = 0; dim < d; dim++) Vs += NE[dim];

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOL_VeInfo = VOLUME->VeInfo;

		ELEMENT = get_ELEMENT_type(VOLUME->type);
		Nve = ELEMENT->Nve;
		for (ve = 0; ve < Nve; ve++) {
			Indve = EToVe[(Vs+(VOLUME->indexg))*NVEMAX+ve];

			for (i = 0; i < NVEINFO; i++)
				VOL_VeInfo[ve+Nve*i] = VeInfo[Indve+NVe*i];
		}
	}
}

static void mark_curved_VOLUME(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d      = DB.d,
	             NVe    = DB.NVe,
	             *NE    = DB.NE,
	             *EToVe = DB.EToVe,
	             *VeInfo = DB.VeInfo;

	// standard datatypes
	unsigned int dim, ve, f, Nve, Vs, *VeInd, NveCurved, NfCurved, fCurved[NFMAX], *VeSurface;

	struct S_ELEMENT *ELEMENT;

	Vs = 0; for (dim = 0; dim < d; dim++) Vs += NE[dim];

	VeSurface = &VeInfo[2*NVe];

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	Nve = ELEMENT->Nve;

	VeInd = VOLUME->VeInd;
	for (ve = 0; ve < Nve; ve++)
		VeInd[ve] = EToVe[(Vs+(VOLUME->indexg))*NVEMAX+ve];

	if (!VOLUME->curved) {
		// Count the number of "curved" vertices and faces in the VOLUME
		for (f = 0; f < NFMAX; f++)
			fCurved[f] = 0;

		NveCurved = 0;
		for (ve = 0; ve < Nve; ve++) {
			if (VeInfo[VeInd[ve]]) {
				fCurved[VeSurface[VeInd[ve]]] = 1;
				NveCurved++;
			}
		}

		NfCurved = 0;
		for (f = 0; f < NFMAX; f++) {
			if (fCurved[f])
				NfCurved++;
		}

		if (NveCurved > NfCurved) {
			if (d != DMAX)
				printf("Error: Should not be entering.\n"), EXIT_MSG;

			VOLUME->curved = 2;
		}
	}

	// Remove vertex index information for straight VOLUMEs
	if (!VOLUME->curved) {
		for (ve = 0; ve < Nve; ve++)
			VeInd[ve] = UINT_MAX;
	}
}

void setup_geometry(void)
{
	// Initialize DB Parameters
	char         *MeshType = DB.MeshType;
	unsigned int ExactGeom = DB.ExactGeom,
	             d         = DB.d,
	             *NE       = DB.NE,
	             *EToVe    = DB.EToVe;
	double       *VeXYZ    = DB.VeXYZ;

	unsigned int PrintTesting = 0;

	// Standard datatypes
	unsigned int ve, dim, P, vn, Vs,
	             NvnGs, NvnGc, NCols;
	double       *XYZ_vC, *XYZ_S, *I_vGs_vGc;

	struct S_ELEMENT   *ELEMENT;
	struct S_VOLUME    *VOLUME;
	struct S_FACET     *FACET;

	// silence
	NvnGs = NvnGc = 0;
	XYZ_S = I_vGs_vGc = NULL;

	Vs = 0; for (dim = 0; dim < d; dim++) Vs += NE[dim];

	// Modify vertex locations if exact geometry is known
	DB.VeInfo = NULL;
	if (ExactGeom) {
		mark_curved_vertices();

		if(!DB.MPIrank && !DB.Testing)
			printf("    Modify vertex nodes if exact geometry is known\n");
		vertices_to_exact_geom();
		initialize_VOLUME_VeInfo();
	}

	// Set up XYZ_S

	/* For the vectorized version of the code, set up a linked list looping through each type of element and each
	 * polynomial order; do this for both VOLUMEs and FACETs (eventually). Then this list can be used to sequentially
	 * performed vectorized operations on each type of element at once.
	 */

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		P      = VOLUME->P;
		XYZ_vC = VOLUME->XYZ_vC;

		ELEMENT = get_ELEMENT_type(VOLUME->type);
		NvnGs = ELEMENT->NvnGs[1];

		// Fix XYZ_vC if modifications were made in projection to exact geometry
		if (ExactGeom && VOLUME->curved) {
			for (ve = 0; ve < NvnGs; ve++) {
			for (dim = 0; dim < d; dim++) {
				XYZ_vC[dim*NvnGs+ve] = VeXYZ[EToVe[(Vs+(VOLUME->indexg))*NVEMAX+ve]*d+dim];
			}}
		}

		if (!VOLUME->curved) {
			// If not curved, the P1 geometry representation suffices to fully specify the element geometry.
			VOLUME->NvnG = NvnGs;

			XYZ_S = malloc(NvnGs*d * sizeof *XYZ_S); // keep
			VOLUME->XYZ_S = XYZ_S;

			for (dim = 0; dim < d; dim++) {
			for (vn = 0; vn < NvnGs; vn++) {
				XYZ_S[dim*NvnGs+vn] = XYZ_vC[dim*NvnGs+vn];
			}}
		} else {
			NvnGc = ELEMENT->NvnGc[P];
			I_vGs_vGc = ELEMENT->I_vGs_vGc[1][P][0];

			NCols = d*1; // d coordinates * 1 element

			VOLUME->NvnG = NvnGc;

			XYZ_S = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep

			mm_d(CBCM,CBT,CBNT,NvnGc,NCols,NvnGs,1.0,0.0,I_vGs_vGc,XYZ_vC,XYZ_S);
		}
		VOLUME->XYZ_S = XYZ_S;
	}

	if (PrintTesting)
		output_to_paraview("ZTest_Geom_straight"); // Output straight coordinates to paraview
//output_to_paraview("ZTest_Geom_straight"); // Output straight coordinates to paraview

	// Set up curved geometry nodes
	if (strstr(MeshType,"ToBeCurved")) {
		if (!DB.MPIrank && !DB.Testing)
			printf("    Set geometry of VOLUME nodes in ToBeCurved Mesh\n");

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			setup_ToBeCurved(VOLUME);
	} else if (strstr(MeshType,"Curved")) {
		// Ensure that VOLUMEs are marked as curved if they only have an edge on a curved boundary
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			mark_curved_VOLUME(VOLUME);
		}

//		if (!DB.MPIrank && !DB.Testing)
		if (!DB.MPIrank)
			printf("    Set geometry of VOLUME nodes in Curved Mesh\n");

		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			setup_Curved(VOLUME);
	} else {
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			setup_straight(VOLUME);
	}

	if (!DB.MPIrank && !DB.Testing)
		printf("    Set FACET XYZ\n");
	for (FACET = DB.FACET; FACET; FACET = FACET->next)
		setup_FACET_XYZ(FACET);

	if (!DB.MPIrank && !DB.Testing)
		printf("    Set up geometric factors\n");
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
		setup_geom_factors(VOLUME);

	if (!DB.MPIrank && !DB.Testing)
		printf("    Set up normals\n");
	for (FACET = DB.FACET; FACET; FACET = FACET->next)
		setup_normals(FACET);

	if (PrintTesting) {
		output_to_paraview("ZTest_Geom_curved"); // Output curved coordinates to paraview
		output_to_paraview("ZTest_Normals");     // Output normals to paraview
	}
}

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const struct S_FACET *FACET,
                     const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int PV, PF, Vtype, Vcurved, FtypeInt;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	PV       = VOLUME->P;
	PF       = FACET->P;
	Vtype    = VOLUME->type;
	Vcurved  = VOLUME->curved;
	FtypeInt = FACET->typeInt;

	ELEMENT     = get_ELEMENT_type(Vtype);
	ELEMENT_OPS = ELEMENT;

	OPS->NvnG = VOLUME->NvnG;
	OPS->NfnS = ELEMENT_OPS->NfnS[PF][IndClass];
	if (FtypeInt == 's') {
		OPS->NfnI = ELEMENT_OPS->NfnIs[PF][IndClass];
		if (!Vcurved) {
			OPS->I_vG_fI = ELEMENT_OPS->I_vGs_fIs[1][PF];
			OPS->I_vG_fS = ELEMENT_OPS->I_vGs_fS[1][PF];
		} else {
			OPS->I_vG_fI = ELEMENT_OPS->I_vGc_fIs[PV][PF];
			OPS->I_vG_fS = ELEMENT_OPS->I_vGc_fS[PV][PF];
		}
	} else {
		OPS->NfnI = ELEMENT_OPS->NfnIc[PF][IndClass];
		if (!Vcurved) {
			OPS->I_vG_fI = ELEMENT_OPS->I_vGs_fIc[1][PF];
			OPS->I_vG_fS = ELEMENT_OPS->I_vGs_fS[1][PF];
		} else {
			OPS->I_vG_fI = ELEMENT_OPS->I_vGc_fIc[PV][PF];
			OPS->I_vG_fS = ELEMENT_OPS->I_vGc_fS[PV][PF];
		}
	}
}
