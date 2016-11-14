// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_Curved.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_FACET.h"

#include "setup_geometry.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "array_norm.h"
#include "array_free.h"

#include "array_print.h" // ToBeDeleted

/*
 *	Purpose:
 *		Set up to be curved meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *		VOLUMEs with curved FACETs have curved = 1 while VOLUMEs with curved EDGEs only have curved = 2.
 *		The operations performed are redundant in 3D. As the blending of FACET to VOLUME is not affected by EDGE node
 *		displacement, EDGE blending is performed for all curved VOLUMEs and the EDGE displacement is recomputed (but not
 *		used) for the FACET node displacement. Check if this is significant while profiling and potentially make
 *		modifications. (ToBeModified)
 *		Potentially add support for different surface parametrizations: Chord Method, Radial Projection, Equal Tangent.
 *
 *	Notation:
 *		XYZ      : High-order VOLUME (XYZ) coordinates after curving.
 *		XYZE_CmS : (XYZ) (E)DGE  coordinates, (C)urved (m)inus (S)traight.
 *		XYZF_CmS : (XYZ) (F)ACET coordinates, (C)urved (m)inus (S)traight.
 *
 *	References:
 */

#define DSPHERE  1
#define CYLINDER 2

struct S_pc {
	unsigned int Nn, VeSurface;
	double       **VeXYZ, *PComps;
};

struct S_XYZ {
	unsigned int Nn, VeSurface;
	double       *XYZ, *PComps;
};

typedef void (*compute_pc_tdef) (struct S_pc *data);
typedef void (*compute_XYZ_tdef) (struct S_XYZ *data);

static void compute_pc_dsphere(struct S_pc *data) {
	/*
	 *	Purpose:
	 *		Compute (t)heta and (p)hi components associated with the input vertices.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn;
	double       r, t, p, *PComps, **VeXYZ;

	if (data->VeSurface == 0)
		r = DB.rIn;
	else if (data->VeSurface == 1)
		r = DB.rOut;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	if (d == 3) {
		printf("Potential modifications needed so that theta is handled correctly at the poles.\n");
	}

	Nn     = data->Nn;
	VeXYZ  = data->VeXYZ;
	PComps = data->PComps;

	for (n = 0; n < Nn; n++) {
		t = atan2(VeXYZ[n][1],VeXYZ[n][0]);

		if (d == 2)
			p = PI/2.0;
		else
			p = acos(VeXYZ[n][2]/r);

		PComps[0*Nn+n] = t;
		PComps[1*Nn+n] = p;
	}
//array_print_d(Nn,2,PComps,'C');
}

static void compute_XYZ_dsphere(struct S_XYZ *data)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn;
	double       *theta, *phi, *XYZ, *X, *Y, *Z, t, p, r;

	if (data->VeSurface == 0)
		r = DB.rIn;
	else if (data->VeSurface == 1)
		r = DB.rOut;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	Nn  = data->Nn;
	XYZ = data->XYZ;

	theta = &data->PComps[0*Nn];
	phi   = &data->PComps[1*Nn];

	X = &XYZ[Nn*0];
	Y = &XYZ[Nn*1];
	if (d == 3)
		Z = &XYZ[Nn*2];

	for (n = 0; n < Nn; n++) {
		t = theta[n];
		p = phi[n];

		X[n] = r*cos(t)*sin(p);
		Y[n] = r*sin(t)*sin(p);
		if (d == 3)
			Z[n] = r*cos(p);
	}
}

static void select_functions_Curved(compute_pc_tdef *compute_pc, compute_XYZ_tdef *compute_XYZ, const unsigned int type)
{
	switch(type) {
	case DSPHERE:
		*compute_pc  = compute_pc_dsphere;
		*compute_XYZ = compute_XYZ_dsphere;
		break;
	default:
		printf("Error: Unsupported.\n"), EXIT_MSG;
		break;
	}
}

void setup_Curved(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *TestCase       = DB.TestCase;
	unsigned int d               = DB.d,
	             *PGc            = DB.PGc,
	             Parametrization = DB.Parametrization;

	// Standard datatypes
	unsigned int dim, e, f, n, ve, Ve, Vf, PV, NvnG, Ne, Nf, Nve, Indve, FuncType,
	             *VeEcon, *VeFcon, Neve, *Nfve, NenG, NfnG, *VeInfo, *VeCurved, *VeSurface, IndSurf,
	             Ecurved, Fcurved, Vcurved;
	double       *XYZ, *XYZ_S, **VeXYZ, *BCoords_vGc, *BCoords_fGc, *BCoords_fGc_dP1, *BCoords_fGs,
	             *PComps_V, *PComps_F, **XYZE_CmS, **XYZF_CmS, *BlendV,
	             *I_vGc_eGc, *I_fGc_vGc, *I_vGc_fGc, *XYZ_update, *XYZ_vC,
	             BlendNum, BlendDen;

	struct S_ELEMENT *ELEMENT, *ELEMENT_F, *ELEMENT_E;
	struct S_pc      *data_pc;
	struct S_XYZ     *data_XYZ;

	// Function pointers
	compute_pc_tdef  compute_pc;
	compute_XYZ_tdef compute_XYZ;

	// Linear portion
	setup_straight(VOLUME);

	Vcurved = VOLUME->curved;
	PV      = VOLUME->P;

	if (!Vcurved || PGc[PV] <= 1)
		return;

	data_pc  = malloc(sizeof *data_pc);  // free
	data_XYZ = malloc(sizeof *data_XYZ); // free

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	NvnG   = VOLUME->NvnG;
	XYZ_S  = VOLUME->XYZ_S;
	XYZ    = VOLUME->XYZ;
	XYZ_vC = VOLUME->XYZ_vC;

	if (strstr(TestCase,"Poisson")) {
		FuncType = DSPHERE;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
	select_functions_Curved(&compute_pc,&compute_XYZ,FuncType);
/*
printf("\n\nVOLUME %d\n",VOLUME->indexg);
array_print_d(NvnG,d,XYZ,'C');
*/
	// Loop over curved EDGEs not on curved FACETs (3D only)
	VeInfo = VOLUME->VeInfo;
	if (d == 3) {
		Ne = ELEMENT->Ne;

		XYZE_CmS = malloc(Ne * sizeof *XYZE_CmS); // free
		for (e = 0; e < Ne; e++) {
			ELEMENT_E = get_ELEMENT_type(LINE);
			NenG = ELEMENT_E->NvnGc[PV];

			XYZE_CmS[e] = calloc(NenG*d , sizeof **XYZE_CmS); // free
		}

		VeEcon = ELEMENT->VeEcon;
		Ne     = ELEMENT->Ne;
		Nve    = ELEMENT->Nve;
		Neve   = NEVEMAX;

		BCoords_vGc = ELEMENT->I_vGs_vGc[1][PV][0];

		VeCurved  = &VeInfo[0*Nve];
		VeSurface = &VeInfo[2*Nve];

		for (e = 0; e < Ne; e++) {
			Ecurved = 1;
			for (ve = 0; ve < Neve; ve++) {
				Indve = VeEcon[e*NEVEMAX+ve];
				if (!ve)
					IndSurf = VeSurface[Indve];

				if (!VeCurved[Indve] || IndSurf != VeSurface[Indve]) {
					Ecurved = 0;
					break;
				}
			}

			if (!Ecurved)
				continue;
printf("sc,e: %d\n",e);
			Ve = e*NEREFMAX;

			ELEMENT_E = get_ELEMENT_type(LINE);
			NenG = ELEMENT_E->NvnGc[PV];

			I_vGc_eGc = ELEMENT->I_vGc_eGc[PV][PV][Ve];

			if (Parametrization == ARC_LENGTH) {
				// Find coordinates of vertices on the EDGE
				VeXYZ = malloc(Neve * sizeof *VeXYZ); // free
				for (ve = 0; ve < Neve; ve++) {
					VeXYZ[ve] = malloc(d * sizeof **VeXYZ); // free
					for (dim = 0; dim < d; dim++)
						VeXYZ[ve][dim] = XYZ_vC[Nve*dim+VeEcon[e*NEVEMAX+ve]];
				}

				// Find barycentric coordinates of VOLUME geometry nodes relating to this EDGE
				BCoords_fGc_dP1 = mm_Alloc_d(CBRM,CBNT,CBNT,NenG,Nve,NvnG,1.0,I_vGc_eGc,BCoords_vGc); // free

				BCoords_fGc = malloc(NenG*Neve * sizeof *BCoords_fGc); // free
				for (n = 0; n < NenG; n++) {
				for (ve = 0; ve < Neve; ve++) {
					BCoords_fGc[n*Neve+ve] = BCoords_fGc_dP1[n*Nve+VeEcon[e*NEVEMAX+ve]];
				}}
				free(BCoords_fGc_dP1);

				// Find values of parametrization components on the EDGE
				PComps_V = malloc(Neve*2 * sizeof *PComps_V); // free

				data_pc->Nn = Neve;
				data_pc->VeSurface = VeInfo[2*Nve+VeEcon[e*NEVEMAX]];
				data_pc->VeXYZ = VeXYZ;
				data_pc->PComps = PComps_V;

				compute_pc(data_pc);
				array_free2_d(Neve,VeXYZ);

				PComps_F = malloc(NenG*2 * sizeof *PComps_F); // free

				mm_CTN_d(NenG,2,Neve,BCoords_fGc,PComps_V,PComps_F);
array_print_d(NenG,2,PComps_F,'C');
EXIT_MSG;
				free(PComps_V);
				free(BCoords_fGc);

				// Compute perturbation of FACET geometry node positions
				data_XYZ->Nn = NenG;
				data_XYZ->VeSurface = VeInfo[2*Nve+VeEcon[e*NEVEMAX]];
				data_XYZ->XYZ = XYZF_CmS[f];
				data_XYZ->PComps = PComps_F;

				compute_XYZ(data_XYZ);
				free(PComps_F);
			}
		}
		array_free2_d(Ne,XYZE_CmS);
	}

	// Loop over curved FACETs
	if (Vcurved == 1) {
		Nf = ELEMENT->Nf;

		XYZF_CmS = malloc(Nf * sizeof *XYZF_CmS); // free
		for (f = 0; f < Nf; f++) {
			ELEMENT_F = get_ELEMENT_F_type(VOLUME->type,f);
			NfnG = ELEMENT_F->NvnGc[PV];

			XYZF_CmS[f] = calloc(NfnG*d , sizeof **XYZF_CmS); // free
		}

		VeFcon      = ELEMENT->VeFcon;
		Nfve        = ELEMENT->Nfve;
		Nve         = ELEMENT->Nve;
		BCoords_vGc = ELEMENT->I_vGs_vGc[1][PV][0];
//array_print_ui(Nve,NVEINFO,VeInfo,'C');

		VeCurved  = &VeInfo[0*Nve];
		VeSurface = &VeInfo[2*Nve];

		for (f = 0; f < Nf; f++) {
			Fcurved = 1;
			for (ve = 0; ve < Nfve[f]; ve++) {
				Indve = VeFcon[f*NFVEMAX+ve];
				if (!ve)
					IndSurf = VeSurface[Indve];

				if (!VeCurved[Indve] || IndSurf != VeSurface[Indve]) {
					Fcurved = 0;
					break;
				}
			}

			if (!Fcurved)
				continue;
//printf("sc: %d\n",f);

			Vf = f*NFREFMAX;

			ELEMENT_F = get_ELEMENT_F_type(VOLUME->type,f);
			NfnG = ELEMENT_F->NvnGc[PV];

			I_fGc_vGc   = ELEMENT->I_fGc_vGc[PV][PV][Vf];
			I_vGc_fGc   = ELEMENT->I_vGc_fGc[PV][PV][Vf];
			BCoords_fGs = ELEMENT->I_fGs_vGc[1][PV][Vf];

// Put this into a separate function later. Figure out what needs to be returned first. ToBeDeleted
			if (Parametrization == ARC_LENGTH) {
				// Find coordinates of vertices on the FACET
				VeXYZ = malloc(Nfve[f] * sizeof *VeXYZ); // free
				for (ve = 0; ve < Nfve[f]; ve++) {
					VeXYZ[ve] = malloc(d * sizeof **VeXYZ); // free
					for (dim = 0; dim < d; dim++)
						VeXYZ[ve][dim] = XYZ_vC[Nve*dim+VeFcon[f*NFVEMAX+ve]];
				}

				// Find barycentric coordinates of VOLUME geometry nodes relating to this FACET
				BCoords_fGc_dP1 = mm_Alloc_d(CBRM,CBNT,CBNT,NfnG,Nve,NvnG,1.0,I_vGc_fGc,BCoords_vGc); // free

				BCoords_fGc = malloc(NfnG*Nfve[f] * sizeof *BCoords_fGc); // free
				for (n = 0; n < NfnG; n++) {
				for (ve = 0; ve < Nfve[f]; ve++) {
					BCoords_fGc[n*Nfve[f]+ve] = BCoords_fGc_dP1[n*Nve+VeFcon[f*NFVEMAX+ve]];
				}}
				free(BCoords_fGc_dP1);

				// Find values of parametrization components on the FACET
				PComps_V = malloc(Nfve[f]*2 * sizeof *PComps_V); // free

				data_pc->Nn = Nfve[f];
				data_pc->VeSurface = VeInfo[2*Nve+VeFcon[f*NFVEMAX]];
				data_pc->VeXYZ = VeXYZ;
				data_pc->PComps = PComps_V;

				compute_pc(data_pc);
				array_free2_d(Nfve[f],VeXYZ);

//array_print_d(NvnG,Nve,BCoords_vGc,'R');
//array_print_d(NfnG,Nfve[f],BCoords_fGc,'R');

				PComps_F = malloc(NfnG*2 * sizeof *PComps_F); // free
				mm_CTN_d(NfnG,2,Nfve[f],BCoords_fGc,PComps_V,PComps_F);
//array_print_d(NfnG,2,PComps_F,'C');
				free(PComps_V);
				free(BCoords_fGc);

				// Compute perturbation of FACET geometry node positions
				data_XYZ->Nn = NfnG;
				data_XYZ->VeSurface = VeInfo[2*Nve+VeFcon[f*NFVEMAX]];
				data_XYZ->XYZ = XYZF_CmS[f];
				data_XYZ->PComps = PComps_F;

				compute_XYZ(data_XYZ);
				free(PComps_F);

				mm_d(CBCM,CBT,CBNT,NfnG,d,NvnG,-1.0,1.0,I_vGc_fGc,XYZ_S,XYZF_CmS[f]);
//array_print_d(NfnG,d,XYZF_CmS[f],'C');

// Make into separate function. ToBeDeleted
				// Compute blend scaling
				// Szabo-Babuska(1991)
				BlendV = malloc(NvnG * sizeof *BlendV); // free

				for (n = 0; n < NvnG; n++) {
					BlendNum = 1.0;
					BlendDen = 1.0;
					for (ve = 0; ve < Nfve[f]; ve++) {
						BlendNum *= BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
						BlendDen *= BCoords_fGs[n*Nfve[f]+ve];
					}
					if (BlendNum < EPS)
						BlendV[n] = 0.0;
					else
						BlendV[n] = BlendNum/BlendDen;
//printf("%d % .3e % .3e\n",n,BlendNum,BlendDen);
				}
//array_print_d(NvnG,1,BlendV,'R');
//array_print_d(NvnG,Nve,BCoords_vGc,'R');
//EXIT_MSG;

//array_print_d(NvnG,Nve,BCoords_vGc,'R');
/*
unsigned int NBI = 15;
double BlendSum, BlendProd;
double *BlendInfo = calloc(NBI*NvnG , sizeof *BlendInfo);
for (n = 0; n < NvnG; n++) {
	double L1, L2, *BlendInfoPtr = &BlendInfo[n*NBI];
	double *I_fGs_vGc = &ELEMENT->I_fGs_vGc[1][PV][Vf][n*2];
	BlendProd = 1.0;
	BlendSum  = 0.0;
	for (ve = 0; ve < Nfve[f]; ve++) {
		BlendProd *= BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
		BlendSum  += BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
	}
	L1 = BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+0]];
	L2 = BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+1]];

	BlendInfoPtr[13] = BlendV[n];
	if (BlendProd < EPS)
		BlendV[n] = 0.0;
	else {
		BlendV[n] = BlendProd/((1.0+2.0*L1*L2-(L1*L1+L2*L2))/4.0);

	BlendInfoPtr[0] = L1;
	BlendInfoPtr[1] = L2;
	BlendInfoPtr[2] = BlendSum;
	BlendInfoPtr[3] = BlendProd;
	BlendInfoPtr[5] = L1*L2/pow((L1+L2),2.0);
	BlendInfoPtr[6] = (1.0+2.0*L1*L2-(L1*L1+L2*L2))/4.0;
	BlendInfoPtr[8] = BlendProd/BlendInfoPtr[5];
	BlendInfoPtr[9] = BlendProd/BlendInfoPtr[6];
	BlendInfoPtr[11] = I_fGs_vGc[0]*I_fGs_vGc[1];
	BlendInfoPtr[11] = L1/(1-L2);
	BlendInfoPtr[14] = BlendV[n];

	BlendV[n] = sqrt(L1/(1-L2)*L2/(1-L1));
	BlendV[n] = (L1/(1-L2)+L2/(1-L1))*0.5;
	double a, b;
	a = L1/(1-L2);
	b = L2/(1-L1);
	BlendV[n] = 2.0*a*b/(a+b);
	BlendV[n] = 2*L1*L2/(L1+L2-(L1*L1+L2*L2));
	BlendV[n] = 4*L1*L2/(1+2*L1*L2-(L1*L1+L2*L2));
	BlendV[n] = BlendSum;
	BlendInfoPtr[11] = BlendV[n];
	}
}
//array_print_d(NvnG,Nve,BCoords_vGc,'R');
array_print_d(NvnG,NBI,BlendInfo,'R');
EXIT_MSG;
*/
				// Blend FACET perturbation to VOLUME geometry nodes
				XYZ_update = mm_Alloc_d(CBCM,CBT,CBNT,NvnG,d,NfnG,1.0,I_fGc_vGc,XYZF_CmS[f]); // free
				for (n = 0; n < NvnG; n++) {
					if (BlendV[n] < EPS)
						continue;

					for (dim = 0; dim < d; dim++)
						XYZ[n+NvnG*dim] += BlendV[n]*XYZ_update[n+NvnG*dim];
				}
				free(XYZ_update);
				free(BlendV);
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}


		}
		array_free2_d(Nf,XYZF_CmS);
	}
/*
array_print_d(NvnG,d,XYZ,'C');
if (Vcurved == 1)
	EXIT_MSG;
*/

	free(data_pc);
	free(data_XYZ);
}
