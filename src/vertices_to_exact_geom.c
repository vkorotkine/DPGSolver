// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "setup_Curved.h"
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

double f_gaussian_bump      (const double x, const double y, const unsigned int d);
double f_naca_symmetric     (const double x, const double y, const unsigned int d);
double f_ellipsoidal_bump   (const double x, const double y, const unsigned int d);
double f_ellipsoidal_corner (const double x, const double y, const unsigned int d);
double f_joukowski_symmetric(const double x, const double y, const unsigned int d);

void select_functions_surface(surface_tdef *f_surface)
{
	// Initialize DB Parameters
	char *Geometry = DB.Geometry;

	if (strstr(Geometry,"GaussianBump")) {
		*f_surface = f_gaussian_bump;
	} else if (strstr(Geometry,"NacaSymmetric")) {
		*f_surface = f_naca_symmetric;
	} else if (strstr(Geometry,"EllipsoidalBump")) {
		*f_surface = f_ellipsoidal_bump;
	} else if (strstr(Geometry,"ExpansionCorner")) {
		*f_surface = f_ellipsoidal_corner;
	} else if (strstr(Geometry,"JoukowskiSymmetric")) {
		*f_surface = f_joukowski_symmetric;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

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
		for (i = 0, iMax = 25; i < iMax; i++) {
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

			if (fabs(f/dfdq) < 3.0*EPS || OnCorner == iMax/2)
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
	if (strstr(Geometry,"GaussianBump")  ||
	    strstr(Geometry,"NacaSymmetric") ||
	    strstr(Geometry,"JoukowskiSymmetric") ||
	    strstr(Geometry,"EllipsoidalBump")) {
		double F_xy;

		surface_tdef f_surface;

		select_functions_surface(&f_surface);

		for (ve = 0; ve < NVe; ve++) {
			if (!VeUpdate[ve])
				continue;

			VeUpdate[ve]  = 0;
			VeSurface[ve] = 0;

			F_xy = f_surface(VeXYZ[ve*d],VeXYZ[ve*d+1],d);
			if (fabs(VeXYZ[ve*d+dM1]-F_xy) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				VeXYZ[ve*d+dM1] = F_xy;
			}
		}
	} else if (strstr(Geometry,"n-Cylinder_Hollow")) {
		unsigned int dCheck = 2;
		double       rIn, rOut, ve_norm2, t;

		rIn  = DB.rIn;
		rOut = DB.rOut;

		for (ve = 0; ve < NVe; ve++) {
			ve_norm2 = array_norm_d(dCheck,&VeXYZ[ve*d],"L2");

			if (fabs(ve_norm2-rIn) < NODETOL_MESH) {
				VeSurface[ve] = 0;
				t = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d]);
				VeXYZ[ve*d]   = rIn*cos(t);
				VeXYZ[ve*d+1] = rIn*sin(t);
			} else if (fabs(ve_norm2-rOut) < NODETOL_MESH) {
				VeSurface[ve] = 1;
				t = atan2(VeXYZ[ve*d+1],VeXYZ[ve*d]);
				VeXYZ[ve*d]   = rOut*cos(t);
				VeXYZ[ve*d+1] = rOut*sin(t);
			}
		}
	} else if (strstr(Geometry,"n-Ball")) {
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
	} else if (strstr(Geometry,"n-Ellipsoid") ||
	           strstr(Geometry,"ExpansionCorner")) {
		double t, p, *abc, X, Y, Z, a, b, c, aIn, aOut;

		aIn  = DB.aIn;
		aOut = DB.aOut;

		abc = malloc(DMAX * sizeof *abc); // free

		for (ve = 0; ve < NVe; ve++) {
			if (!VeUpdate[ve])
				continue;

			VeUpdate[ve] = 0;

			X = VeXYZ[ve*d+0];
			Y = VeXYZ[ve*d+1];
			Z = VeXYZ[ve*d+(d-1)];

			get_abc_ellipse(1,(double *) &VeXYZ[ve*d],abc);

			a = abc[0];
			b = abc[1];
			c = abc[2];

			if (fabs(a-aIn) < EPS) {
				VeSurface[ve] = 0;
			} else if (fabs(a-aOut) < EPS) {
				VeSurface[ve] = 1;
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}

			t = atan2(Y/b,X/a);

			if (d == 2) {
				p = PI/2.0;
			} else {
				if (Z-c > 1.0)
					p = acos(1.0);
				else
					p = acos(Z/c);
			}

			VeXYZ[ve*d+0] = a*cos(t)*sin(p);
			VeXYZ[ve*d+1] = b*sin(t)*sin(p);
			if (d == 3)
				VeXYZ[ve*d+2] = c*cos(p);
		}
		free(abc);
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
	} else if (strstr(Geometry,"n-Cube")) {
		// Do nothing
	} else if (strstr(Geometry,"HoldenRamp")) {
		if (!strstr(Geometry,"HoldenRampCurved")) {
			// Do nothing
		} else {
			// Initialize DB Parameters
			double rIn = DB.rIn;

			// Standard datatypes
			double Xc, Yc, XcR, t, x, y;

			rIn = DB.rIn;

			Xc = -rIn/tan(82.5/180.0*PI);
			Yc =  rIn;
			XcR = Xc + rIn*cos(3.0/2.0*PI+15.0/180.0*PI);

			for (ve = 0; ve < NVe; ve++) {
				if (!VeUpdate[ve])
					continue;

				VeUpdate[ve]  = 0;
				VeSurface[ve] = 0;

				// Correct y-coordinate based on which surface the vertex is located.
				x = VeXYZ[ve*d];
				y = VeXYZ[ve*d+1];

				// Treat straight line segments
				if (x-Xc < EPS) {
//					printf("ls: % .3e % .3e % .3e\n",x,Xc,VeXYZ[ve*d+1]);
					VeXYZ[ve*d+1] = 0.0;
					continue;
				} else if (x-XcR > EPS) {
//					printf("rs: % .3e % .3e % .3e\n",x,XcR,VeXYZ[ve*d+1]-tan(15.0/180.0*PI)*VeXYZ[ve*d+0]);
					VeXYZ[ve*d+1] = tan(15.0/180.0*PI)*VeXYZ[ve*d+0];
					continue;
				}

				t = atan2(y-Yc,x-Xc);
				while (t < 0.0)
					t += 2*PI;
				if (t < 3.0/2.0*PI || t > 3.0/2.0*PI+15.0/180.0*PI)
					printf("Error: Invalid value (t = % .3e).\n",t), EXIT_MSG;

//printf("%d % .13e % .13e % .13e\n",ve,t,t-1.5*PI,t-3.0/2.0*PI+15.0/180.0*PI);

				// Radial projection to the boundary
				VeXYZ[ve*d+0] = Xc+rIn*cos(t);
				VeXYZ[ve*d+1] = Yc+rIn*sin(t);
				// No need to correct z component
			}
		}
	} else {
		printf("%s\n",Geometry);
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

double f_gaussian_bump(const double x, const double y, const unsigned int d)
{
	// ToBeModified: Note that this function is replicated in EvalTPFunction
	double a, b, c;

	if (d == 2) {
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

double f_naca_symmetric(const double x, const double y, const unsigned int d)
{
	unsigned int i;
	double c, t, a[5], Output;

	c = DB.NSc;
	t = DB.NSt;

	a[0] = DB.NS0;
	a[1] = DB.NS1;
	a[2] = DB.NS2;
	a[3] = DB.NS3;
	a[4] = DB.NS4;

	Output = a[0]*sqrt(x/c);
	for (i = 1; i < 5; i++)
		Output += a[i]*pow(x/c,(double) i);
	Output *= 5.0*c*t;

	if (x < EPS || x > 1.0-EPS)
		Output = 0.0;

	return Output;

	if (0) // silence
		printf("%d %f\n",d,y);
}

double f_ellipsoidal_bump(const double x, const double y, const unsigned int d)
{
	double a, b, c;

	if (d == 2) {
		a = DB.aIn;
		b = DB.bIn;

		if (fabs(x) > a-EPS)
			return 0.0;
		else
			return b*sqrt(1.0-pow(x/a,2.0));
	} else if (d == 3) {
		c = DB.cIn;
		printf("Add support.\n"), EXIT_MSG;
		printf("%e %e\n",y,c); // silence
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

double f_ellipsoidal_corner(const double x, const double y, const unsigned int d)
{
	if (d == 2) {
		return f_ellipsoidal_bump(x,y,d);
	} else if (d == 3) {
		// Ensure that the same function can be used in 3D as well if this becomes supported.
		printf("Add support.\n"), EXIT_MSG;
		printf("%e\n",y); // silence
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

double f_joukowski_symmetric(const double x, const double y, const unsigned int d)
{
	// Initialize DB Parameters
	double a  = DB.JSa,
	       l  = DB.JSl,
	       xL = DB.JSxL;

	// Standard datatypes
	double t, p, l2, l3;

	l2 = pow(l,2.0);
	l3 = pow(l,3.0);

	// Find parametrization coordinates ((t)heta, (p)hi)
	if (d == 2) {
		unsigned int count, countMax;
		double       cost, cost2, f, dfdt;

		if (x < xL+EPS || x > 2*a-EPS) {
//			printf("Return0 (% .3e % .3e)\n",x,xL);
			return 0.0;
		}

		// Use Newton's method to find t corresponding to the input x
		t = 0.5*PI;
		countMax = 100;
		for (count = 0; count < countMax; count++) {
			cost  = cos(t);
			cost2 = pow(cost,2.0);
			f    = -2.0*a*(l3+(l3+2.0*l2+l)*cost2+l2-(2.0*l3+3.0*l2+2.0*l+1.0)*cost+l)/
			       (2.0*l2-2.0*(l2+l)*cost+2.0*l+1.0)-x;
			dfdt = 4.0*(l3+(l3+2.0*l2+l)*cost2+l2-(2.0*l3+3.0*l2+2.0*l+1.0)*cost+l)*(l2+l)*sin(t)/
			       pow(2.0*l2-2.0*(l2+l)*cost+2.0*l+1.0,2.0)
			     + 2.0*(2.0*(l3+2.0*l2+l)*cost*sin(t)-(2.0*l3+3.0*l2+2.0*l+1.0)*sin(t))/
			       (2.0*l2-2.0*(l2+l)*cost+2.0*l+1.0);

			if (fabs(dfdt) > EPS)
				t -= f/dfdt;
			else
				t = PI;

			if (t < 0.0)
				t = 0.0;
			else if (t > PI)
				t = PI;

			if (fabs(f) < 1e1*EPS)
				break;
		}

		if (count == countMax) {
			// Attempt to find the point using the bisection method
			unsigned int Found;
			double       dt, sign_t, t2, f2;

			if (t < 0.5*PI) {
				dt = t;
				sign_t = -1.0;
			} else {
				dt = PI-t;
				sign_t = 1.0;
			}

			Found = 0;
			t2 = t;
			for (count = 0; count < 50; count++) {
				cost  = cos(t2+sign_t*dt);
				cost2 = pow(cost,2.0);

				f2 = -2.0*a*(l3+(l3+2.0*l2+l)*cost2+l2-(2.0*l3+3.0*l2+2.0*l+1.0)*cost+l)/
			        (2.0*l2-2.0*(l2+l)*cost+2.0*l+1.0)-x;

				if (sign(f) == sign(f2))
					t2 += sign_t*dt;

				if (fabs(f2) < 1e1*EPS) {
					Found = 1;
					break;
				}

				dt /= 2;
			}

			if (Found) {
				t = t2;
			} else {
				// If still not found
				printf("% .3e % .3e % .3e\n",t,t2,sign_t);
				printf("Error: Newton's Method not converging (% .3e % .3e % .3e % .3e % .3e)\n",x,t,t-PI,f,dfdt),
			    EXIT_MSG;
			}
		}

		if (t < SQRT_EPS || t > PI-SQRT_EPS)
			return 0.0;
		else
			return 2.0*a*(l3+2.0*l2-(l3+2.0*l2+l)*cos(t)+l)*sin(t)/(2.0*l2-2.0*(l2+l)*cos(t)+2.0*l+1.0);
	} else if (d == 3) {
		printf("Add support.\n"), EXIT_MSG;
		p = 0.0;
		printf("%f %f\n",y,p);
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}
