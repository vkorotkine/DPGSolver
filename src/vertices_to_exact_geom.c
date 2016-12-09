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
	char         *Geometry = DB.Geometry;
	double       *VeXYZ    = DB.VeXYZ;

	// Standard datatypes
	unsigned int dM1, ve, *VeUpdate, *VeSurface;

	VeUpdate  = &VeInfo[1*NVe];
	VeSurface = &VeInfo[2*NVe];

// Removed case for TestCase == Test (which did nothing). Ensure that vertices are not moved when they should not be for
// the test runs. Or just check that all tests are passing. (ToBeDeleted)

	dM1 = d-1;
	if (strstr(Geometry,"GaussianBump")) {
		double F_xy;

		for (ve = 0; ve < NVe; ve++) {
			F_xy = f_gaussian_bump(VeXYZ[ve*d],VeXYZ[ve*d+1],d);
			if (fabs(VeXYZ[ve*d+dM1]-F_xy) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				VeXYZ[ve*d+dM1] = F_xy;
			}
		}
	} else if (strstr(Geometry,"Annular_Section")) {
		double rIn, rOut, ve_norm2, t;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");

			if (fabs(ve_norm2-rIn) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				t = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(t);
				VeXYZ[ve*d+1] = rIn*sin(t);
			} else if (fabs(ve_norm2-rOut) < NODETOL_MESH) {
				VeSurface[ve] = 1;
				t = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d+2]);
				VeXYZ[ve*d]   = rIn*cos(t);
				VeXYZ[ve*d+1] = rOut*sin(t);
			}
		}
	} else if (strstr(Geometry,"dm1-Spherical_Section")) {
		double rIn, rOut, r, ve_norm2, t, p;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < NVe; ve++) {
			if (!VeUpdate[ve])
				continue;

			ve_norm2 = array_norm_d(d,&VeXYZ[ve*d],"L2");

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
	} else if (strstr(Geometry,"Ringleb")) {
		unsigned int *VeRingleb;
		double       Q0, KMin, KMax, q, k, x, y, a, rho, J;

		Q0   = DB.Q0;
		KMin = DB.KMin;
		KMax = DB.KMax;

		VeRingleb = &VeInfo[3*NVe];

		for (ve = 0; ve < NVe; ve++) {
			if (!VeUpdate[ve])
				continue;

			VeUpdate[ve] = 0;

			// Find which surface this vertex is located on.
			if (VeRingleb[ve] == 'f') {
				x = VeXYZ[ve*d];
				y = VeXYZ[ve*d+1];

				q = Q0;
				if (y < 0.0) {
					VeSurface[ve] = 0; // Inflow
				} else {
					VeSurface[ve] = 1; // Outflow
				}

				a   = sqrt(1.0-0.5*GM1*q*q);
				rho = pow(a,2.0/GM1);
				J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

				k = sqrt(2.0/(2.0*rho*(x+0.5*J)+1.0/(q*q)));

				// Correct y-coordinate
				y = sign(y)/(k*rho*q)*sqrt(1.0-pow(q/k,2.0));
				VeXYZ[ve*d+1] = y;
			} else if (VeRingleb[ve] == 'w') {
				unsigned int count, countMax;
				double       sign_y, f, dfdq, relax;

				x = VeXYZ[ve*d];
				y = VeXYZ[ve*d+1];

				if (y < 0.0)
					sign_y = -1.0;
				else
					sign_y =  1.0;

				if (x < 0.0) {
					k = KMax;
					VeSurface[ve] = 2; // Inner wall
				} else {
					k = KMin;
					VeSurface[ve] = 3; // Outer wall
				}

				// Use Newton's method (with relaxation) to find q (using equation for y-coordinate)
				q = 0.5*(Q0+k);

				countMax = 50;
				relax = 0.5;
				for (count = 0; count < countMax; count++) {
					a   = sqrt(1.0-0.5*GM1*q*q);
					rho = pow(a,2.0/GM1);

					f = sign_y/(k*k*rho*q)*sqrt(k*k-q*q)-y;
					// Do not compute denominator if f == 0.0 as it may be 'inf'
					if (fabs(f) < EPS)
						dfdq = 1.0;
					else
						dfdq = sign_y*(pow((1.0-0.5*GM1*q*q),-1.0/GM1)*(2.0*pow(q,4.0)-k*k*((GAMMA+1)*q*q-2.0))/
						               (k*k*q*q*(GM1*q*q-2.0)*sqrt(k*k-q*q)));
					q -= relax*f/dfdq;

					if (q < Q0)
						q = Q0;
					else if (q > k)
						q = k;

					relax = (1-fabs(f/dfdq))*1.0;

					if (fabs(f/dfdq) < 1e2*EPS)
						break;
				}
				if (count == countMax)
					printf("Error: Newton's method not converging.\n"), EXIT_MSG;

				a   = sqrt(1.0-0.5*GM1*q*q);
				rho = pow(a,2.0/GM1);
				J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

				// Correct x-coordinate
				x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-0.5*J;
				VeXYZ[ve*d] = x;
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}
		}
EXIT_MSG;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void vertices_to_exact_geom_VOLUME(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *Geometry = DB.Geometry;
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

	if (strstr(Geometry,"dm1-Spherical_Section")) {
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

		a = DB.GBa;
		b = DB.GBb;
		c = DB.GBc;

		return a*exp(-pow(x-b,2)/(2.0*pow(c,2)));
	} else if (d == 3) {
		printf("Set up the 3D Gaussian bump case.\n"), exit(1);
		// some function of x and y
		printf("%e\n",y);
	} else {
		printf("Error: Invalid dimension used in f_gaussian_bump.\n"), exit(1);
	}
}
