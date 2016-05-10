#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Move vertices to correspond with the exact geometry (eliminating errors from the mesh generator)
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
*/

double f_gaussian_bump(const double x, const double y, const unsigned int d);

void vertices_to_exact_geom(void)
{
	// Initialize DB Parameters
	unsigned int d         = DB.d,
	             NVe       = DB.NVe;
	char         *TestCase = DB.TestCase;
	double       *VeXYZ    = DB.VeXYZ;

	// Standard datatypes
	unsigned int dM1, ve;

	dM1 = d-1;
	if (strstr(TestCase,"GaussianBump") != NULL) {
		double F_xy;

		for (ve = 0; ve < NVe; ve++) {
			F_xy = f_gaussian_bump(VeXYZ[ve*d],VeXYZ[ve*d+1],d);
			if (fabs(VeXYZ[ve*d+dM1]-F_xy) < NODETOL_MESH)
				VeXYZ[ve*d+dM1] = F_xy;
		}
	} else if (strstr(TestCase,"SupersonicVortex") != NULL) {
		double rIn, rOut, ve_norm2, theta;

		rIn  = 1.0;
		rOut = 1.384;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");

			if (fabs(ve_norm2-rIn) < NODETOL_MESH) {
				theta = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(theta);
				VeXYZ[ve*d+1] = rIn*sin(theta);
			} else if (fabs(ve_norm2-rOut) < NODETOL_MESH) {
				theta = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(theta);
				VeXYZ[ve*d+1] = rOut*sin(theta);
			}
		}
	} else if (strstr(TestCase,"dSphericalBump") != NULL) {
		double r, ve_norm2, theta, phi;

		r = 0.1;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");
			if (fabs(ve_norm2-r) < NODETOL_MESH) {
				theta = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d]);

				VeXYZ[ve*d]   = r*cos(theta);
				VeXYZ[ve*d+1] = r*sin(theta);

				if (d == 3) {
					phi = acos(VeXYZ[ve*d+2]/r);

					VeXYZ[ve*d]   *= sin(phi);
					VeXYZ[ve*d+1] *= sin(phi);
					VeXYZ[ve*d+2]  = r*cos(phi);
				}
			}
		}
	}
}

double f_gaussian_bump(const double x, const double y, const unsigned int d)
{
	// ToBeModified: Note that this function is replicated in EvalTPFunction
	double a, b, c;

	if (d == 2) {
		printf("Make sure that the same bump parameters are being used across all functions. Exiting.\n"), exit(1);
		// use y when d = 2 to eliminate compiler warning
		a = y;

		a = 0.0625;
		b = 0.0;
		c = 0.2;

		return a*exp(-pow(x-b,2)/(2.0*pow(c,2)));
	} else if (d == 3) {
		printf("Set up the 3D Gaussian bump case.\n"), exit(1);
		// some function of x and y
	} else {
		printf("Error: Invalid dimension used in f_gaussian_bump.\n"), exit(1);
	}
}
