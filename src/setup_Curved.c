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

#include "element_functions.h"
#include "matrix_functions.h"
#include "array_norm.h"

#include "array_print.h" // ToBeDeleted

/*
 *	Purpose:
 *		Set up to be curved meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *		VOLUMEs with curved FACETs have curved = 1 while VOLUMEs with curved EDGEs only have curved = 2.
 *		Potentially add support for different surface parametrizations: Chord Method, Radial Projection, Equal Tangent.
 *
 *	Notation:
 *		XYZ : High-order VOLUME (XYZ) coordinates after curving.
 *
 *	References:
 */

#define DSPHERE_Poisson 1
#define DSPHERE_Euler   2
#define CYLINDER        3

struct S_pc {
	unsigned int Nn, VeSurface;
	double       **VeXYZ, *PComps, rIn, rOut;
};

struct S_XYZ {
	unsigned int Nn, VeSurface, *Node_on_F;
	double       *XYZ, *PComps, rIn, rOut;
};

typedef void (*compute_pc_tdef) (struct S_pc *data);
typedef void (*compute_XYZ_tdef) (struct S_XYZ *data);

static void compute_pc_dsphere(struct S_pc *data) {
	/*
	 *	Purpose:
	 *		Compute (t)heta and (p)hi components associated with the input vertices.
	 *
	 *	Comments:
	 *		Need to modify this so that theta at the poles of the sphere is handled correctly. (ToBeDeleted)
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn;
	double       r, t, p, *PComps, **VeXYZ;

	if (data->VeSurface == 0)
		r = data->rIn;
	else if (data->VeSurface == 1)
		r = data->rOut;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	if (d == 3) {
		printf("See comments.\n");
		printf("Error: Unsupported.\n"), EXIT_MSG;
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
}

static void compute_XYZ_dsphere(struct S_XYZ *data)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn, *Node_on_F;
	double       *theta, *phi, *XYZ, *X, *Y, *Z, t, p, r;

	if (data->VeSurface == 0)
		r = data->rIn;
	else if (data->VeSurface == 1)
		r = data->rOut;
	else
		printf("Error: Unsupported.\n"), EXIT_MSG;

	Nn        = data->Nn;
	Node_on_F = data->Node_on_F;
	XYZ       = data->XYZ;

	theta = &data->PComps[0*Nn];
	phi   = &data->PComps[1*Nn];

	X = &XYZ[Nn*0];
	Y = &XYZ[Nn*1];
	if (d == 3)
		Z = &XYZ[Nn*2];

	for (n = 0; n < Nn; n++) {
		if (!Node_on_F[n])
			continue;

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
	case DSPHERE_Poisson:
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
	             NVe             = DB.NVe,
	             Parametrization = DB.Parametrization,
	             *VeInfo         = DB.VeInfo;
	double       *VeXYZ          = DB.VeXYZ;

	// Standard datatypes
	unsigned int i, f, n, ve, PV, NvnG, Vcurved, Nf, Nve, *VeFcon, *Nfve, *VeInd, *Node_on_F, FuncType, *Fmask;
	double       *XYZ, *XYZ_S, **VeFXYZ, *BCoords_vGc, *BCoords_fGc, *PComps_V, *PComps_F;

	struct S_ELEMENT *ELEMENT;
	struct S_FACET   *FACET;
	struct S_pc      *data_pc;
	struct S_XYZ     *data_XYZ;

	// Function pointers
	compute_pc_tdef  compute_pc;
	compute_XYZ_tdef compute_XYZ;

	data_pc  = malloc(sizeof *data_pc);  // free
	data_XYZ = malloc(sizeof *data_XYZ); // free

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	NvnG  = VOLUME->NvnG;
	XYZ_S = VOLUME->XYZ_S;

	XYZ = malloc (NvnG*d * sizeof *XYZ); // keep
	VOLUME->XYZ = XYZ;

	if (strstr(TestCase,"Poisson")) {
		FuncType = DSPHERE_Poisson;
		data_pc->rIn  = data_XYZ->rIn  = 0.5;
		data_pc->rOut = data_XYZ->rOut = 1.0;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
	select_functions_Curved(&compute_pc,&compute_XYZ,FuncType);

	// ToBeDeleted
	for (i = 0; i < NvnG*d; i++)
		XYZ[i] = XYZ_S[i];

	/*
	 *	Strategy:
	 *	1) Loop over EDGEs (3D only):
	 *		- Find all curved EDGEs which are not on curved FACETs.
	 *			Use VeInfo for this.
	 *		- Project straight geometry to surface.
	 *			Arclength Parametrization: Given two components for the surface parametrization (sphere: theta, phi;
	 *			cylinder: theta, z), use the end point values and interpolate between them according to the spacing of
	 *			the EDGE geometry nodes. Use the Barycentric coordinate projection of the parametrization components to
	 *			define the surface points (Also works for FACETs).
	 *				*** Careful of undefined theta at sphere poles. Can deal with this by considering a sphere which is
	 *				rotated about the phi-axis such that no points lie at phi = 0 or PI. ***
	 *			Radial Projection: Simple.
	 *			Chord Method (Trasdahl(2012)): 2D description given here (standard extension to 3D). Find the equation
	 *			of the line connecting two face vertices in physical space. Project straight lines normal to this line
	 *			which touch the exact surface. This gives the values to be used at the GLL nodes.
	 *			Equal Tangent Parametrization (Trasdahl(2012)): Compute the nodes on the surface using Newton's method.
	 *			See the thesis for details.
	 *		- Blend curved EDGE geometry to adjacent FACETs.
	 *	2) Loop over FACETs:
	 *		- Find all curved FACETs.
	 *		- Project straight geometry to surface.
	 *			See description above.
	 *	3) Blend FACET information to the VOLUME.
	 */

	Vcurved = VOLUME->curved;

	// Loop over curved EDGEs not on curved FACETs (3D only)
	if (d == 3 && Vcurved == 2) {
		printf("Error: Add support.\n"), EXIT_MSG;
	}

	// Loop over curved FACETs
	if (Vcurved == 1) {
		PV    = VOLUME->P;
		VeInd = VOLUME->VeInd;
		NvnG  = VOLUME->NvnG;

		Nf          = ELEMENT->Nf;
		VeFcon      = ELEMENT->VeFcon;
		Nfve        = ELEMENT->Nfve;
		Nve         = ELEMENT->Nve;
		BCoords_vGc = ELEMENT->I_vGs_vGc[1][PV][0];

		Node_on_F = malloc(NvnG * sizeof *Node_on_F); // free
		for (f = 0; f < Nf; f++) {
			FACET = VOLUME->FACET[f*NSUBFMAX];

			if (!FACET->curved)
				continue;

// Put this into a separate function later. Figure out what needs to be returned first. ToBeDeleted
			if (Parametrization == ARC_LENGTH) {
				// Find coordinates of vertices on the FACET
				VeFXYZ = malloc(Nfve[f] * sizeof *VeFXYZ); // free
				for (ve = 0; ve < Nfve[f]; ve++)
					VeFXYZ[ve] = &VeXYZ[VeInd[VeFcon[f*NFVEMAX+ve]]*d];

				// Find barycentric coordinates of VOLUME geometry nodes relating to this FACET
				BCoords_fGc = malloc(NvnG*Nfve[f] * sizeof *BCoords_fGc); // free
				for (n = 0; n < NvnG; n++) {
					for (ve = 0; ve < Nfve[f]; ve++) {
						BCoords_fGc[n*Nfve[f]+ve] = BCoords_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
					}

					// Mark nodes which are on the FACET
					if (fabs(array_norm_d(Nfve[f],&BCoords_fGc[n*Nfve[f]],"L1") - 1.0) < EPS)
						Node_on_F[n] = 1;
					else
						Node_on_F[n] = 0;
				}

				// Find values of parametrization components on the FACET
				PComps_V = malloc(Nfve[f]*2 * sizeof *PComps_V); // free

				data_pc->Nn = Nfve[f];
				data_pc->VeSurface = VeInfo[2*NVe+VeInd[VeFcon[f*NFVEMAX]]];
				data_pc->VeXYZ = VeFXYZ;
				data_pc->PComps = PComps_V;

				compute_pc(data_pc);
				free(VeFXYZ);

array_print_d(NvnG,Nve,BCoords_vGc,'R');
array_print_d(NvnG,Nfve[f],BCoords_fGc,'R');

				PComps_F = malloc(NvnG*2 * sizeof *PComps_F); // free
				mm_CTN_d(NvnG,2,Nfve[f],BCoords_fGc,PComps_V,PComps_F);
				free(PComps_V);
				free(BCoords_fGc);

				// Compute new positions for FACET geometry nodes
				data_XYZ->Nn = NvnG;
				data_XYZ->VeSurface = VeInfo[2*NVe+VeInd[VeFcon[f*NFVEMAX]]];
				data_XYZ->XYZ = XYZ;
				data_XYZ->PComps = PComps_F;
				data_XYZ->Node_on_F = Node_on_F;

				compute_XYZ(data_XYZ);
				free(PComps_F);

array_print_ui(NvnG,1,Node_on_F,'R');
			Fmask = ELEMENT->Fmask[PV][PV][f*NFREFMAX];
array_print_ui(PV+1,1,Fmask,'R');
EXIT_MSG;
/*
printf("sC: %d %d %d\n",VOLUME->indexg,FACET->indexg,PV);
for (ve = 0; ve < Nfve[f]; ve++)
	array_print_d(1,d,VeFXYZ[ve],'R');

array_print_d(NvnG,d,XYZ_S,'C');
array_print_d(Nfve[f],2,PComps_V,'C');
array_print_d(NvnG,2,PComps_F,'C');
array_print_d(NvnG,d,XYZ,'C');
exit(1);
*/
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}


		}
		free(Node_on_F);
	}

	free(data_pc);
}
