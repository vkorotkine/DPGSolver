// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "vertices_to_exact_geom.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "element_functions.h"
#include "array_norm.h"

#include "array_print.h"

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
	             NVe       = DB.NVe,
				 *VeInfo   = DB.VeInfo;
	char         *TestCase = DB.TestCase;
	double       *VeXYZ    = DB.VeXYZ;

	// Standard datatypes
	unsigned int dM1, ve, *VeUpdate, *VeSurface;

	VeUpdate = &VeInfo[1*NVe];
	VeSurface= &VeInfo[2*NVe];

	dM1 = d-1;
	if (strstr(TestCase,"GaussianBump")) {
		double F_xy;

		for (ve = 0; ve < NVe; ve++) {
			F_xy = f_gaussian_bump(VeXYZ[ve*d],VeXYZ[ve*d+1],d);
			if (fabs(VeXYZ[ve*d+dM1]-F_xy) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				VeXYZ[ve*d+dM1] = F_xy;
			}
		}
	} else if (strstr(TestCase,"SupersonicVortex")) {
		double rIn, rOut, ve_norm2, theta;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");

			if (fabs(ve_norm2-rIn) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				theta = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(theta);
				VeXYZ[ve*d+1] = rIn*sin(theta);
			} else if (fabs(ve_norm2-rOut) < NODETOL_MESH) {
				VeSurface[ve] = 1;
				theta = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(theta);
				VeXYZ[ve*d+1] = rOut*sin(theta);
			}
		}
	} else if (strstr(TestCase,"dSphericalBump")) {
		double r, ve_norm2, theta, phi;

		r = DB.rIn;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");
			if (fabs(ve_norm2-r) < NODETOL_MESH) {
				VeSurface[ve] = 0;
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
	} else if (strstr(TestCase,"Poisson")) {
		double rIn, rOut, r, ve_norm2, t, p;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");

			if (!VeUpdate[ve])
				continue;

			VeUpdate[ve] = 0;

			if (fabs(ve_norm2-rIn) < NODETOL_MESH) {
				r = rIn;
				VeSurface[ve] = 0;
			} else if (fabs(ve_norm2-rOut) < NODETOL_MESH) {
				r = rOut;
				VeSurface[ve] = 1;
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}

			t = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d]);

			if (d == 2) {
				VeXYZ[ve*d]   = r*cos(t);
				VeXYZ[ve*d+1] = r*sin(t);
			} else if (d == 3) {
				p = acos(VeXYZ[ve*d+2]/r);

				VeXYZ[ve*d]   = r*cos(t)*sin(p);
				VeXYZ[ve*d+1] = r*sin(t)*sin(p);
				VeXYZ[ve*d+2] = r*cos(p);
			}
		}
	} else if (strstr(TestCase,"Test")) {
		// No vertex movement required.
	}
}

void vertices_to_exact_geom_VOLUME(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int ve, Nve, *VeUpdate, *VeSurface;
	double       *XYZ, *X, *Y, *Z;

	struct S_ELEMENT *ELEMENT;

	// silence
	Z = NULL;

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	Nve = ELEMENT->Nve;

	VeUpdate  = &VOLUME->VeInfo[Nve*1];
	VeSurface = &VOLUME->VeInfo[Nve*2];

	XYZ = VOLUME->XYZ_vV;

	X = &XYZ[0*Nve];
	Y = &XYZ[1*Nve];
	if (d == 3)
		Z = &XYZ[2*Nve];

	if (strstr(TestCase,"Poisson")) {
		double rIn, rOut, r, rXYZ, t, p;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < Nve; ve++) {
			if (!VeUpdate[ve])
				continue;

			if (VeSurface[ve] == 0)
				r = rIn;
			else if (VeSurface[ve] == 1)
				r = rOut;
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;

			t = atan2(Y[ve],X[ve]);

			if (d == 2) {
				X[ve] = r*cos(t);
				Y[ve] = r*sin(t);
			} else if (d == 3) {
				rXYZ = sqrt(X[ve]*X[ve]+Y[ve]*Y[ve]+Z[ve]*Z[ve]);
				p = acos(Z[ve]/rXYZ);

				X[ve] = r*cos(t)*sin(p);
				Y[ve] = r*sin(t)*sin(p);
				Z[ve] = r*cos(p);
			}
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
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
