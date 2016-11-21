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
#include "S_FACE.h"

#include "setup_geometry.h"
#include "element_functions.h"
#include "matrix_functions.h"
#include "array_norm.h"
#include "array_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Set up "Curved" meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *		VOLUMEs with curved FACEs have curved = 1 while VOLUMEs with curved EDGEs only have curved = 2.
 *		The operations performed are redundant in 3D. As the blending of FACE to VOLUME is not affected by EDGE node
 *		displacement, EDGE blending is performed for all curved VOLUMEs and the EDGE displacement is recomputed (but not
 *		used) for the FACE node displacement. Check if this is significant while profiling and potentially make
 *		modifications. (ToBeModified)
 *		Potentially add support for different surface parametrizations: Chord Method, Radial Projection, Equal Tangent.
 *
 *	Notation:
 *		XYZ      : High-order VOLUME (XYZ) coordinates after curving.
 *		XYZ_CmS : (XYZ) BOUNDARY (EDGE or FACE) coordinates, (C)urved (m)inus (S)traight.
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

struct S_Blend {
	unsigned int b, NvnG, NbnG, Nve, *Nbve, *VeBcon, NbveMax, *VeInfo, EclassV, type;
	double       *XYZ, *XYZ_vV,
	             *I_vGs_vGc, *I_bGs_vGc, *I_vGc_bGc, *I_bGc_vGc;

	compute_pc_tdef  compute_pc;
	compute_XYZ_tdef compute_XYZ;
};

static void compute_pc_dsphere(struct S_pc *data) {
	/*
	 *	Purpose:
	 *		Compute (t)heta and (p)hi parametrization components associated with the input vertices.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn, OnPole = 0, IndPole;
	double       r, t, p, t_avg, *PComps, **VeXYZ, X, Y, Z, rXYZ;

	if (data->VeSurface == 0)
		r = DB.rIn;
	else if (data->VeSurface == 1)
		r = DB.rOut;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	Nn     = data->Nn;
	VeXYZ  = data->VeXYZ;
	PComps = data->PComps;

	for (n = 0; n < Nn; n++) {
		X = VeXYZ[n][0];
		Y = VeXYZ[n][1];

		t = atan2(Y,X);

		if (d == 2) {
			p = PI/2.0;
		} else {
			Z    = VeXYZ[n][2];
			rXYZ = sqrt(X*X+Y*Y+Z*Z);

			p = acos(Z/rXYZ);
		}

		PComps[0*Nn+n] = t;
		PComps[1*Nn+n] = p;
	}

	if (d == 3) {
		// If one (and only one) of the vertices lies on a pole of the sphere, re-calculate theta for the vertex based
		// on the theta values of the other vertices.
		for (n = 0; n < Nn; n++) {
			p = acos(VeXYZ[n][2]/r);

			// Note: SQRT_EPS used here because cos(x) ~= 1 - 0.5*x^2 for x << 1
			if (fabs(p-0.0) < SQRT_EPS || fabs(p-PI) < SQRT_EPS) {
				OnPole++;
				IndPole = n;
			}
		}

		if (!OnPole) {
			// Do nothing
		} else if (OnPole == 1) {
			t_avg = 0.0;
			for (n = 0; n < Nn; n++) {
				if (n == IndPole)
					continue;
				t      = atan2(VeXYZ[n][1],VeXYZ[n][0]);
				t_avg += t;
			}
			t_avg /= (Nn-1);

			PComps[0*Nn+IndPole] = t_avg;
		} else if (OnPole == 2) {
			printf("Error: Unsupported.\n");
		}
	}
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

static double *compute_BlendV(struct S_Blend *data)
{
	/*
	 *	Purpose:
	 *		Compute the scaling to be used for each of the VOLUME geometry nodes based on the boundary perturbation due
	 *		to the curved geometry.
	 *
	 *	Comments:
	 *		SZABO_BABUSKA blending is based on a generalization of eq. (6.22) in Szabo(1991). Note that the projected
	 *		nodes used along the curved face must be defined as in setup_blending (setup_operators).
	 *
	 *	References:
	 *		Szabo(1991)-Finite_Element_Analysis
	 */

	// Initialize DB Parameters
	unsigned int Blending = DB.Blending;

	// Standard datatypes
	unsigned int b, n, ve, NvnG, Nve, *Nbve, *VeBcon, NbveMax, EclassV, type;
	double       *I_vGs_vGc, *I_bGs_vGc, *BlendV, BlendNum, BlendDen, Blend_yes;

	NvnG = data->NvnG;
	BlendV = malloc(NvnG * sizeof *BlendV); // keep

	type    = data->type;
	EclassV = data->EclassV;

	b       = data->b;
	Nve     = data->Nve;
	Nbve    = data->Nbve;
	VeBcon  = data->VeBcon;
	NbveMax = data->NbveMax;

	I_vGs_vGc = data->I_vGs_vGc;
	I_bGs_vGc = data->I_bGs_vGc;

	if (Blending == HESTHAVEN && type == TRI) {
		for (n = 0; n < NvnG; n++) {
			BlendNum = I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+0]];
			BlendDen = 1.0-I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+1]];
			if (BlendNum < EPS)
				BlendV[n] = 0.0;
			else
				BlendV[n] = BlendNum/BlendDen;
		}
	} else if (Blending == NIELSON && type == TRI) {
		for (n = 0; n < NvnG; n++) {
			// Ensure that no blending is performed for nodes on the FACEs
			Blend_yes = 1;
			for (ve = 0; ve < Nbve[b]; ve++) {
				if (fabs(I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]]) < EPS) {
					Blend_yes = 0;
					break;
				}
			}

			if (Blend_yes)
				BlendV[n] = pow(1.0-I_vGs_vGc[n*Nve+b],2.0);
			else
				BlendV[n] = 0.0;
		}
	} else if (Blending == SZABO_BABUSKA || EclassV == C_SI) {
		for (n = 0; n < NvnG; n++) {
			BlendNum = 1.0;
			BlendDen = 1.0;
			for (ve = 0; ve < Nbve[b]; ve++) {
				BlendNum *= I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]];
				BlendDen *= I_bGs_vGc[n*Nbve[b]+ve];
			}
			if (BlendNum < EPS)
				BlendV[n] = 0.0;
			else
				BlendV[n] = BlendNum/BlendDen;
		}
	} else if (Blending == GORDON_HALL && EclassV == C_TP) {
		for (n = 0; n < NvnG; n++) {
			BlendV[n] = 0.0;
			for (ve = 0; ve < Nbve[b]; ve++)
				BlendV[n] += I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]];
		}
	} else {
		printf("Error: Unsupported.\n");
	}

if (b == 2) {
printf("%d\n",b);
array_print_d(NvnG,1,BlendV,'R');
array_print_d(NvnG,Nve,I_vGs_vGc,'R');
array_print_d(NvnG,Nbve[b],I_bGs_vGc,'R');

		for (n = 0; n < NvnG; n++) {
			BlendNum = 1.0;
			BlendDen = 1.0;
			for (ve = 0; ve < Nbve[b]; ve++) {
				BlendNum *= I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]];
				BlendDen *= I_bGs_vGc[n*Nbve[b]+ve];
			}
			if (BlendNum < EPS)
				BlendV[n] = 0.0;
			else
				BlendV[n] = BlendNum/BlendDen;
		}


array_print_d(NvnG,1,BlendV,'R');
EXIT_MSG;
}

	return BlendV;
}

static double *compute_XYZ_update(struct S_Blend *data)
{
	/*
	 *	Purpose:
	 *		Compute update to VOLUME XYZ geometry components based on the selected curved surface parametrization.
	 */

	// Initialize DB Parameters
	unsigned int d               = DB.d,
	             Parametrization = DB.Parametrization;

	// Standard datatypes
	unsigned int b, n, ve, dim, NvnG, NbnG, Nve, *Nbve, *VeBcon, *VeInfo, NbveMax;
	double       **VeXYZ, *XYZ, *XYZ_vV, *PComps_V, *PComps_B, *XYZ_CmS, *XYZ_update, *XYZ_S,
	             *I_vGs_bGc, *I_vGc_bGc, *I_vGs_vGc, *I_bGs_bGc, *I_bGc_vGc;

	struct S_pc  *data_pc;
	struct S_XYZ *data_XYZ;

	// Function pointers
	compute_pc_tdef  compute_pc;
	compute_XYZ_tdef compute_XYZ;

	NbnG = data->NbnG;
	XYZ_CmS = calloc(NbnG*d , sizeof *XYZ_CmS); // free

	NvnG    = data->NvnG;
	b       = data->b;
	Nve     = data->Nve;
	VeInfo  = data->VeInfo;
	VeBcon  = data->VeBcon;
	NbveMax = data->NbveMax;

	XYZ       = data->XYZ;
	I_vGc_bGc = data->I_vGc_bGc;
	I_bGc_vGc = data->I_bGc_vGc;

	data_pc  = malloc(sizeof *data_pc);  // free
	data_XYZ = malloc(sizeof *data_XYZ); // free

	compute_pc  = data->compute_pc;
	compute_XYZ = data->compute_XYZ;

	if (Parametrization == ARC_LENGTH) {
		Nbve    = data->Nbve;

		XYZ_vV    = data->XYZ_vV;
		I_vGs_vGc = data->I_vGs_vGc;

		// Find coordinates of vertices on the BOUNDARY
		VeXYZ = malloc(Nbve[b] * sizeof *VeXYZ); // free
		for (ve = 0; ve < Nbve[b]; ve++) {
			VeXYZ[ve] = malloc(d * sizeof **VeXYZ); // free
			for (dim = 0; dim < d; dim++)
				VeXYZ[ve][dim] = XYZ_vV[Nve*dim+VeBcon[b*NbveMax+ve]];
		}

		// Find barycentric coordinates of VOLUME geometry nodes relating to this BOUNDARY
		I_vGs_bGc = mm_Alloc_d(CBRM,CBNT,CBNT,NbnG,Nve,NvnG,1.0,I_vGc_bGc,I_vGs_vGc); // free

		I_bGs_bGc = malloc(NbnG*Nbve[b] * sizeof *I_bGs_bGc); // free
		for (n = 0; n < NbnG; n++) {
		for (ve = 0; ve < Nbve[b]; ve++) {
			I_bGs_bGc[n*Nbve[b]+ve] = I_vGs_bGc[n*Nve+VeBcon[b*NbveMax+ve]];
		}}
		free(I_vGs_bGc);

		// Find values of parametrization components on the BOUNDARY
		PComps_V = malloc(Nbve[b]*2 * sizeof *PComps_V); // free

		data_pc->Nn        = Nbve[b];
		data_pc->VeXYZ     = VeXYZ;
		data_pc->PComps    = PComps_V;
		data_pc->VeSurface = VeInfo[2*Nve+VeBcon[b*NbveMax]];

		compute_pc(data_pc);
		array_free2_d(Nbve[b],VeXYZ);

		PComps_B = malloc(NbnG*2 * sizeof *PComps_B); // free
		mm_CTN_d(NbnG,2,Nbve[b],I_bGs_bGc,PComps_V,PComps_B);
		free(PComps_V);
		free(I_bGs_bGc);

		// Compute perturbation of BOUNDARY geometry node positions
		data_XYZ->Nn        = NbnG;
		data_XYZ->XYZ       = XYZ_CmS;
		data_XYZ->PComps    = PComps_B;
		data_XYZ->VeSurface = VeInfo[2*Nve+VeBcon[b*NbveMax]];

		compute_XYZ(data_XYZ);
		free(PComps_B);

		mm_d(CBCM,CBT,CBNT,NbnG,d,NvnG,-1.0,1.0,I_vGc_bGc,XYZ,XYZ_CmS);
	} else if (Parametrization == RADIAL_PROJECTION) {
		// Compute XYZ_C
		XYZ_S = malloc(NbnG*d * sizeof *XYZ_S); // free

		mm_d(CBCM,CBT,CBNT,NbnG,d,NvnG,1.0,0.0,I_vGc_bGc,XYZ,XYZ_S);

		VeXYZ = malloc(NbnG   * sizeof *VeXYZ); // free
		for (n = 0; n < NbnG; n++) {
			VeXYZ[n] = malloc(d * sizeof *VeXYZ[n]); // free
			for (dim = 0; dim < d; dim++)
				VeXYZ[n][dim] = XYZ_S[n+NbnG*dim];
		}

		PComps_B = malloc(NbnG*2 * sizeof *PComps_B); // free
		
		data_pc->Nn        = NbnG;
		data_pc->VeXYZ     = VeXYZ;
		data_pc->PComps    = PComps_B;
		data_pc->VeSurface = VeInfo[2*Nve+VeBcon[b*NbveMax]];

		compute_pc(data_pc);
		array_free2_d(NbnG,VeXYZ);

		data_XYZ->Nn        = NbnG;
		data_XYZ->XYZ       = XYZ_CmS;
		data_XYZ->PComps    = PComps_B;
		data_XYZ->VeSurface = VeInfo[2*Nve+VeBcon[b*NbveMax]];

		compute_XYZ(data_XYZ);
		free(PComps_B);

		// Subtract XYZ_S;
		for (n = 0; n < NbnG*d; n++)
			XYZ_CmS[n] -= XYZ_S[n];

		free(XYZ_S);
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	free(data_pc);
	free(data_XYZ);

	XYZ_update = mm_Alloc_d(CBCM,CBT,CBNT,NvnG,d,NbnG,1.0,I_bGc_vGc,XYZ_CmS); // free
	free(XYZ_CmS);

	return XYZ_update;
}

static void blend_boundary(struct S_VOLUME *VOLUME, const unsigned int BType)
{
	/*
	 *	Purpose:
	 *		Compute difference between straight geometry and that of the parametrized curved surface and blend it into
	 *		the VOLUME, also adding additional Vertex contributions if applicable.
	 *
	 *	Comments:
	 *		The (b)oundary prefix is used here in place of either (e)dge or (f)ace, depending on BType.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             Blending  = DB.Blending;

	// Standard datatypes
	unsigned int b, ve, dim, n, Vb, PV, Indve, IndSurf, Vtype,
	             Nve, Nb, *Nbve, *VeBcon, NbveMax, NbrefMax, NvnG, Bcurved, FuncType,
	             *VeInfo, *VeCurved, *VeSurface;
	double       *XYZ, *XYZ_update, *XYZ_vV, *BlendV, *BlendVe, *I_vGs_vGc;

	struct S_ELEMENT *ELEMENT, *ELEMENT_B;
	struct S_pc      *data_pc;
	struct S_XYZ     *data_XYZ;
	struct S_Blend   *data_blend;

	// silence
	Nb = NbrefMax = NbveMax = IndSurf = 0;
	VeBcon = Nbve = NULL;
	ELEMENT_B = NULL;

	data_pc    = malloc(sizeof *data_pc);    // free
	data_XYZ   = malloc(sizeof *data_XYZ);   // free
	data_blend = malloc(sizeof *data_blend); // free

	if (strstr(TestCase,"Poisson")) {
		FuncType = DSPHERE;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
	select_functions_Curved(&data_blend->compute_pc,&data_blend->compute_XYZ,FuncType);

	PV     = VOLUME->P;
	Vtype  = VOLUME->type;
	NvnG   = VOLUME->NvnG;
	XYZ    = VOLUME->XYZ;
	XYZ_vV = VOLUME->XYZ_vV;
	VeInfo = VOLUME->VeInfo;

	data_blend->NvnG   = NvnG;
	data_blend->XYZ    = XYZ;
	data_blend->XYZ_vV = XYZ_vV;
	data_blend->VeInfo = VeInfo;

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	data_blend->type = Vtype;

	Nve       = ELEMENT->Nve;
	I_vGs_vGc = ELEMENT->I_vGs_vGc[1][PV][0];

	data_blend->EclassV   = ELEMENT->Eclass;
	data_blend->Nve       = Nve;
	data_blend->I_vGs_vGc = I_vGs_vGc;

	VeCurved  = &VeInfo[0*Nve];
	VeSurface = &VeInfo[2*Nve];

	if (BType == 'e') {
		Nb     = ELEMENT->Ne;
		Nbve   = malloc(Nb * sizeof *Nbve); // free
		for (b = 0; b < Nb; b++)
			Nbve[b] = ELEMENT->Neve;
		VeBcon = ELEMENT->VeEcon;

		NbveMax  = NEVEMAX;
		NbrefMax = NEREFMAX;
	} else if (BType == 'f') {
		Nb     = ELEMENT->Nf;
		Nbve   = ELEMENT->Nfve;
		VeBcon = ELEMENT->VeFcon;

		NbveMax  = NFVEMAX;
		NbrefMax = NFREFMAX;
	}

	data_blend->Nbve    = Nbve;
	data_blend->NbveMax = NbveMax;
	data_blend->VeBcon  = VeBcon;

	for (b = 0; b < Nb; b++) {
		Bcurved = 1;
		for (ve = 0; ve < Nbve[b]; ve++) {
			Indve = VeBcon[b*NbveMax+ve];
			if (!ve)
				IndSurf = VeSurface[Indve];

			if (!VeCurved[Indve] || IndSurf != VeSurface[Indve]) {
				Bcurved = 0;
				break;
			}
		}

		if (!Bcurved)
			continue;

		Vb = b*NbrefMax;
		if (BType == 'e') {
			ELEMENT_B = get_ELEMENT_type(LINE);
			data_blend->I_vGc_bGc = ELEMENT->I_vGc_eGc[PV][PV][Vb];
			data_blend->I_bGc_vGc = ELEMENT->I_eGc_vGc[PV][PV][Vb];
			data_blend->I_bGs_vGc = ELEMENT->I_eGs_vGc[1][PV][Vb];
		} else if (BType == 'f') {
			ELEMENT_B = get_ELEMENT_F_type(ELEMENT->type,b);
			data_blend->I_vGc_bGc = ELEMENT->I_vGc_fGc[PV][PV][Vb];
			data_blend->I_bGc_vGc = ELEMENT->I_fGc_vGc[PV][PV][Vb];
			data_blend->I_bGs_vGc = ELEMENT->I_fGs_vGc[1][PV][Vb];
		}

		data_blend->b    = b;
		data_blend->NbnG = ELEMENT_B->NvnGc[PV];

		BlendV     = compute_BlendV(data_blend);     // free
		XYZ_update = compute_XYZ_update(data_blend); // free

		// Blend BOUNDARY perturbation to VOLUME geometry nodes
		for (n = 0; n < NvnG; n++) {
			if (BlendV[n] < EPS)
				continue;

			for (dim = 0; dim < d; dim++)
				XYZ[n+NvnG*dim] += BlendV[n]*XYZ_update[n+NvnG*dim];
		}
		free(XYZ_update);
		free(BlendV);
	}

	// Add additional contribution from vertices if applicable
	if (Blending == NIELSON && Vtype == TRI) {
		BlendVe = malloc(NvnG*Nve * sizeof *BlendVe); // free

		for (n = 0; n < NvnG; n++) {
		for (ve = 0; ve < Nve; ve++) {
			BlendVe[n*Nve+ve] = pow(I_vGs_vGc[n*Nve+ve],2.0);
		}}
		mm_d(CBCM,CBT,CBNT,NvnG,d,Nve,0.0,1.0,BlendVe,XYZ_vV,XYZ);
		free(BlendVe);
	} else {
		// Do nothing
	}

	free(data_pc);
	free(data_XYZ);
	free(data_blend);

	if (BType == 'e') {
		free(Nbve);
	}
}

void setup_Curved(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             *PGc = DB.PGc;

	// Standard datatypes
	unsigned int Vcurved;

	// Linear portion
	setup_straight(VOLUME);

	Vcurved = VOLUME->curved;
	if (!Vcurved || PGc[VOLUME->P] <= 1)
		return;

	// Treat curved EDGEs not on curved FACEs (3D only)
	if (d == DMAX)
		blend_boundary(VOLUME,'e');

	// Treat curved FACEs
	if (Vcurved == 1)
		blend_boundary(VOLUME,'f');
}
