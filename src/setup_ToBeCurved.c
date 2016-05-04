#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Set up to be curved meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *
 *	Notation:
 *
 *	References:
 *		Gordon(1973)_Construction of Curvilinear Coordinate Systems and Applications to Mesh Generation
 *		Gordon(1973)_Transfinite Element Methods- Blending-Function Interpolation over Arbitraty Curved Element Domains
 *		Rosca(2011)_Uniform Spherical Grids via Equal Area Projection from the Cube to the Sphere
 */

void setup_ToBeCurved(void)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, dim,
	             NvnG;
	double *XYZ, *XYZs;

	struct S_VOLUME *VOLUME;

	if (strstr(TestCase,"dSphericalBump")   != NULL ||
		strstr(TestCase,"PorousdSphere")    != NULL ||
		strstr(TestCase,"SupersonicVortex") != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			NvnG = VOLUME->NvnG;
			XYZs = VOLUME->XYZs;

			XYZ = malloc (NvnG*d * sizeof *XYZ); // keep

			VOLUME->XYZ = XYZ;
		}
	} else if (strstr(TestCase,"GaussianBump")     != NULL ||
	           strstr(TestCase,"PolynomialBump")   != NULL) {
		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			NvnG = VOLUME->NvnG;
			XYZs = VOLUME->XYZs;

			XYZ = malloc (NvnG*d * sizeof *XYZ); // keep

			VOLUME->XYZ = XYZ;
		}
	} else if (strstr(TestCase,"PeriodicVortex") != NULL) {
		double n = 2.0, A = 0.1, L0 = 2.0, dxyz = 1.0;
		double *X0, *Y0, *Z0;

		// silence
		Z0 = NULL;

		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			NvnG = VOLUME->NvnG;
			XYZs = VOLUME->XYZs;

			XYZ = malloc (NvnG*d * sizeof *XYZ); // keep

			// d > 1 for this case
			for (dim = 0; dim < d; dim++) {
				X0 = &XYZs[0*NvnG];
				Y0 = &XYZs[1*NvnG];
				if (dim == 2)
					Z0 = &XYZs[2*NvnG];
			}

			if (d == 2) {
				for (i = 0; i < NvnG; i++) {
					XYZ[       i] = X0[i] + A*dxyz*sin(n*PI/L0*Y0[i]);
					XYZ[1*NvnG+i] = Y0[i] + A*dxyz*sin(n*PI/L0*X0[i]);
				}
			} else if (d == 3) {
				for (i = 0; i < NvnG; i++) {
					XYZ[       i] = X0[i] + A*dxyz*sin(n*PI/L0*Y0[i])*sin(n*PI/L0*Z0[i]);
					XYZ[1*NvnG+i] = Y0[i] + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*Z0[i]);
					XYZ[2*NvnG+i] = Z0[i] + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*Y0[i]);
				}
			} else {
				printf("Error: PeriodicVortex TestCase not supported for dimension d = %d.\n",d), exit(1);
			}
			VOLUME->XYZ = XYZ;
		}
	} else {
		printf("Error: Unsupported TestCase for the ToBeCurved MeshType.\n"), exit(1);
	}
}
