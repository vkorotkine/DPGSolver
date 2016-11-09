// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "setup_Curved.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "S_DB.h"
#include "S_VOLUME.h"
/*
#include <math.h>

#include "Parameters.h"
#include "Macros.h"

#include "array_norm.h"
#include "cubature.h"

#include "array_print.h" // ToBeDeleted
*/

/*
 *	Purpose:
 *		Set up to be curved meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *
 *	Notation:
 *		XYZ : High-order VOLUME (XYZ) coordinates after curving.
 *
 *	References:
 */

void setup_Curved(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, NvnG;
	double       *XYZ, *XYZ_S;

	NvnG  = VOLUME->NvnG;
	XYZ_S = VOLUME->XYZ_S;

	XYZ = malloc (NvnG*d * sizeof *XYZ); // keep
	VOLUME->XYZ = XYZ;

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
	 *		- Blend curved EDGE geometry to adjacent FACETs.
	 *	2) Loop over FACETs:
	 *		- Find all curved FACETs.
	 *		- Project straight geometry to surface.
	 */

	if (strstr(TestCase,"Poisson")) {
//		rIn  = 0.5;
//		rOut = 1.0;
	}
}
