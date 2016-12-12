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
#include "matrix_functions.h"
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

void Ringleb_boundary(double *xStore, double *yStore, double qIn, double kIn, const char RinglebType)
{
	unsigned int OnCorner = 0;
	double       a, rho, J, x, y, sign_y;

	x = *xStore;

	sign_y = 1.0;
	if (*yStore < 0.0)
		sign_y = -1.0;

	if (RinglebType == 'f') {
		const double q = qIn;
		double       k;

		a   = sqrt(1.0-0.5*GM1*q*q);
		rho = pow(a,2.0/GM1);
		J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

		k = sqrt(2.0/(2.0*rho*(x+0.5*J)+1.0/(q*q)));
		y = sign_y/(k*k*rho*q)*sqrt(k*k-q*q);
	} else if (RinglebType == 'w') {
		// Initialize DB Parameters
		double Q0 = DB.Q0;

		// Standard datatypes
		unsigned int i, iMax;
		const double k = kIn;
		double       q, f, dadq, drhodq, dJdq, dfdq;

		// Use Newton's method to find q (using equation for x-coordinate)
		q = 0.5*(k+Q0);
		OnCorner = 0;
		for (i = 0, iMax = 20; i < iMax; i++) {
			a   = sqrt(1.0-0.5*GM1*q*q);
			rho = pow(a,2.0/GM1);
			J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

			f = 0.5/rho*(2.0/(k*k)-1.0/(q*q))-0.5*J-x;

			dadq   = 0.5*pow(1.0-0.5*GM1*q*q,-0.5)*(-GM1*q);
			drhodq = 2.0/GM1*pow(a,2.0/GM1-1.0)*dadq;
			dJdq   = -1.0*(pow(a,-2.0)+pow(a,-4.0)+pow(a,-6.0)-1.0/(a*a-1.0))*dadq;
			dfdq   = -pow(k*rho,-2.0)*drhodq-0.5*(-pow(rho*q*q,-2.0)*(q*q*drhodq+2*rho*q)+dJdq);

			q -= f/dfdq;
			if (q < Q0 || fabs(q-Q0) < EPS) {
				OnCorner++;
				q = Q0;
			} else if (q > k || fabs(q-k) < EPS) {
				OnCorner++;
				q = k;
			}

			if (fabs(f/dfdq) < 1e2*EPS || OnCorner == iMax/2)
				break;
		}
		if (i == iMax)
			printf("Error: Newton's method not converging (% .3e).\n",fabs(f/dfdq)), EXIT_MSG;

		if (!OnCorner) {
			y = sign_y/(k*k*rho*q)*sqrt(k*k-q*q);
		} else {
			// Fix x-coordinate as well
			a   = sqrt(1.0-0.5*GM1*q*q);
			rho = pow(a,2.0/GM1);
			J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

			x = 0.5/rho*(2.0/(k*k)-1.0/(q*q))-0.5*J;
			y = sign_y/(k*k*rho*q)*sqrt(k*k-q*q);
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	// Correct y-coordinate
	*xStore = x;
	*yStore = y;
}

void vertices_to_exact_geom(void)
{
	/*
	 *	Purpose:
	 *		Move vertices from the initial mesh file to ensure that they are on the boundary to within close to machine
	 *		precision.
	 *
	 *	Comments:
	 *		The surface parametrization used to achieve this purpose here is not so crucial as the points are already
	 *		very close to being correctly placed on the boundary. However, using an analogous function to move newly
	 *		created vertices to the boundary must more carefully account for a good surface parametrization so as to
	 *		choose good vertex placement in regions of high curvature.
	 */

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
		double       Q0, KMin, KMax, k, x, y;

		Q0   = DB.Q0;
		KMin = DB.KMin;
		KMax = DB.KMax;

		VeRingleb = &VeInfo[3*NVe];

		for (ve = 0; ve < NVe; ve++) {
			if (!VeUpdate[ve])
				continue;

			VeUpdate[ve] = 0;

			// Correct y-coordinate based on which surface the vertex is located.
			x = VeXYZ[ve*d];
			y = VeXYZ[ve*d+1];

			if (VeRingleb[ve] == 'f') {
				if (y > 0.0) {
					VeSurface[ve] = 0; // Inflow
				} else {
					VeSurface[ve] = 1; // Outflow
				}

				Ringleb_boundary(&VeXYZ[ve*d+0],&VeXYZ[ve*d+1],Q0,0.0,'f');
			} else if (VeRingleb[ve] == 'w') {
				if (x < 0.0) {
					k = KMax;
					VeSurface[ve] = 2; // Inner wall
				} else {
					k = KMin;
					VeSurface[ve] = 3; // Outer wall
				}

				Ringleb_boundary(&VeXYZ[ve*d+0],&VeXYZ[ve*d+1],0.0,k,'w');
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}
		}
	} else if (strstr(Geometry,"PeriodicVortex")) {
		// Do nothing
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
