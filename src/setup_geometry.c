// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

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
 *		XYZ_S : High-order VOLUME (XYZ) coordinates at (S)tart (i.e. before curving)
 *		XYZ   : High-order VOLUME (XYZ) coordinates after curving (if applicable)
 *
 *	References:
*/

void setup_geometry(void)
{
	// Initialize DB Parameters
	char         *MeshType = DB.MeshType,
	             *TestCase = DB.TestCase;
	unsigned int ExactGeom = DB.ExactGeom,
	             d         = DB.d,

	             Testing   = DB.Testing;

	int          PrintTesting = 0;

	// Standard datatypes
	unsigned int dim, P, vn,
	             NvnG, NvnGs, NvnGc, NCols;
	double       *XYZ_vC, *XYZ_S, *XYZ,
	             *I_vGs_vGc;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;
	struct S_FACET   *FACET;

	// silence
	NvnGs = 0; NvnGc = 0;
	XYZ_S = NULL;
	I_vGs_vGc = NULL;

	// Modify vertex locations if exact geometry is known
	if (ExactGeom) {
		if(!DB.MPIrank) printf("    Modify vertex nodes if exact geometry is known\n");
		printf("Did not yet verify the implementation.\n");
		vertices_to_exact_geom();
	}

	// Set up XYZ_S

	/* For the vectorized version of the code, set up a linked list looping through each type of element and each
	 * polynomial order; do this for both VOLUMEs and FACETs (eventually). Then this list can be used to sequentially
	 * performed vectorized operations on each type of element at once.
	 */

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		P      = VOLUME->P;
		XYZ_vC = VOLUME->XYZ_vC;

		ELEMENT = get_ELEMENT_type(VOLUME->type);
		if (!VOLUME->curved) {
			// If not curved, the P1 geometry representation suffices to fully specify the element geometry.
			NvnG = ELEMENT->NvnGs[0];

			VOLUME->NvnG = NvnG;

			XYZ_S = malloc(NvnG*d * sizeof *XYZ_S); // keep
			XYZ   = malloc(NvnG*d * sizeof *XYZ);   // keep
			VOLUME->XYZ_S = XYZ_S;

			for (dim = 0; dim < d; dim++) {
			for (vn = 0; vn < NvnG; vn++) {
				XYZ_S[dim*NvnG+vn] = XYZ_vC[dim*NvnG+vn];
			}}
		} else {
			NvnGs = ELEMENT->NvnGs[0];
			NvnGc = ELEMENT->NvnGc[P];
			I_vGs_vGc = ELEMENT->I_vGs_vGc[P];

			NCols = d*1; // d coordinates * 1 element

			VOLUME->NvnG = NvnGc;

			XYZ_S = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
			XYZ   = malloc(NvnGc*NCols * sizeof *XYZ);   // keep

			mm_d(CblasColMajor,CblasTrans,CblasNoTrans,NvnGc,NCols,NvnGs,1.0,I_vGs_vGc,XYZ_vC,XYZ_S);
		}
		VOLUME->XYZ_S = XYZ_S;

//array_print_d(VOLUME->NvnG,d,VOLUME->XYZ_S,'C');
//exit(1);
	}

	if (Testing)
		output_to_paraview("ZTest_Geom_straight"); // Output straight coordinates to paraview

	// Set up curved geometry nodes
	if (strstr(MeshType,"ToBeCurved") != NULL) {
		printf("    Set geometry of VOLUME nodes in ToBeCurved Mesh\n");
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
			setup_ToBeCurved(VOLUME);
	} else {
		printf("Add in support for MeshType != ToBeCurved");
		exit(1);
	}

	printf("    Set up geometric factors\n");
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		setup_geom_factors(VOLUME);

	printf("    Set up normals\n");
	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next)
		setup_normals(FACET);

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next)
		free(VOLUME->C_vC);

	if (Testing) {
		output_to_paraview("ZTest_Geom_curved"); // Output curved coordinates to paraview
		output_to_paraview("ZTest_Normals");     // Output normals to paraview
	}
}
