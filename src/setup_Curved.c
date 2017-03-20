// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

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
#include "vertices_to_exact_geom.h"
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
 *	References: ToBeModified
 *		Szabo(1991)-Finite_Element_Analysis (eq. 6.22)
 *		Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (section 9.1.1)
 *		Nielson(1977)-The_Side-Vertex_Method_for_Interpolation_in_Triangles (eq. 2.2, 2.13)
 *		Gordon(1973)-Transfinite_Element_Methods-_Blending-Function_Interpolation_over_Arbitrary_Curved_Element_Domains
 */

struct S_pc {
	unsigned int Nn, VeSurface;
	double       **VeXYZ, *PComps;
};

struct S_XYZ {
	unsigned int Nn, VeSurface;
	double       *XYZ, *PComps;
};

typedef void   (*compute_pc_tdef)  (struct S_pc *data);
typedef void   (*compute_XYZ_tdef) (struct S_XYZ *data);

struct S_Blend {
	unsigned int b, NvnG, NbnG, Nve, *Nbve, *VeBcon, NbveMax, *VeInfo, EclassV, type, BC;
	double       *XYZ, *XYZ_vV,
	             *I_vGs_vGc, *I_bGs_vGc, *I_vGc_bGc, *I_bGc_vGc;

	compute_pc_tdef  compute_pc;
	compute_XYZ_tdef compute_XYZ;
};

void compute_plane(const double *XYZ1, const double *XYZ2, const double *XYZ3, double *n, double *d_p)
{
	/*
	 *	Purpose:
	 *		Compute normal vector to a plane defined by three points.
	 *
	 *	Comments:
	 *		The plane is defined by: a*x+b*y+c*z = d, where the n = (a,b,c).
	 */

	unsigned int i, d;
	double       Vec1[3], Vec2[3];

	d = 3;

	for (i = 0; i < d; i++) {
		Vec1[i] = XYZ1[i]-XYZ3[i];
		Vec2[i] = XYZ2[i]-XYZ3[i];
	}

	// compute cross product
	n[0] =  (Vec1[1]*Vec2[2]-Vec1[2]*Vec2[1]);
	n[1] = -(Vec1[0]*Vec2[2]-Vec1[2]*Vec2[0]);
	n[2] =  (Vec1[0]*Vec2[1]-Vec1[1]*Vec2[0]);

	*d_p = 0.0;
	for (i = 0; i < d; i++)
		*d_p += n[i]*XYZ3[i];
}

void get_abc_ellipse(const unsigned int Nn, double *XYZ, double *abc)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d    = DB.d;
	double       rIn  = DB.rIn,
	             aIn  = DB.aIn,
	             bIn  = DB.bIn,
	             cIn  = DB.cIn,
	             rOut = DB.rOut,
	             aOut = DB.aOut,
	             bOut = DB.bOut,
	             cOut = DB.cOut;

	// Standard datatypes
	unsigned int n;
	double       *x, *y, *z, norm_In, norm_Out;

	if (strstr(TestCase,"PrandtlMeyer")) {
		abc[0] = aIn;
		abc[1] = bIn;
		abc[2] = cIn;
		return;
	}

	x = &XYZ[0*Nn];
	y = &XYZ[1*Nn];
	z = &XYZ[(d-1)*Nn];

	norm_In  = 0.0;
	norm_Out = 0.0;

	for (n = 0; n < Nn; n++) {
		if (d == 2) {
			norm_In  += fabs(sqrt(pow(x[n]/aIn,2.0) +pow(y[n]/bIn,2.0))-1.0);
			norm_Out += fabs(sqrt(pow(x[n]/aOut,2.0)+pow(y[n]/bOut,2.0))-1.0);
		} else {
			norm_In  += fabs(sqrt(pow(x[n]/aIn,2.0) +pow(y[n]/bIn,2.0) +pow(z[n]/cIn,2.0)) -1.0);
			norm_Out += fabs(sqrt(pow(x[n]/aOut,2.0)+pow(y[n]/bOut,2.0)+pow(z[n]/cOut,2.0))-1.0);
		}
	}
	norm_In  /= Nn;
	norm_Out /= Nn;

	if (norm_In < 4e-1*(rOut-rIn)) {
		abc[0] = aIn;
		abc[1] = bIn;
		abc[2] = cIn;
	} else if (norm_Out < 4e-1*(rOut-rIn)) {
		abc[0] = aOut;
		abc[1] = bOut;
		abc[2] = cOut;
	} else {
		array_print_d(Nn,d,XYZ,'C');
		printf("% .3e % .3e\n",norm_In,norm_Out);
		printf("Error: Did not find the ellipse.\n"), EXIT_MSG;
	}
}

static double get_radius(const unsigned int Nn, double *XYZ)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d;
	double       rIn  = DB.rIn,
	             rOut = DB.rOut;

	// Standard datatypes
	unsigned int n;
	double       *x, *y, *z, norm_rIn, norm_rOut, r, r2;

	x = &XYZ[0*Nn];
	y = &XYZ[1*Nn];
	z = &XYZ[(d-1)*Nn];

	norm_rIn  = 0.0;
	norm_rOut = 0.0;

	for (n = 0; n < Nn; n++) {
		r2 = x[n]*x[n]+y[n]*y[n];
		if (d == 3)
			r2 += z[n]*z[n];
		r          = sqrt(r2);
		norm_rIn  += sqrt((r-rIn)* (r-rIn));
		norm_rOut += sqrt((r-rOut)*(r-rOut));
	}
	norm_rIn  /= Nn;
	norm_rOut /= Nn;

	if (norm_rIn < 4e-1*(rOut-rIn)) {
		r = rIn;
	} else if (norm_rOut < 4e-1*(rOut-rIn)) {
		r = rOut;
	} else {
		printf("% .3e % .3e\n",norm_rIn,norm_rOut);
		array_print_d(Nn,d,XYZ,'C');
		printf("Error: Did not find the radius.\n"), EXIT_MSG;
	}
	return r;
}

void compute_normal_displacement(const unsigned int Nn, const unsigned int curved_normal, const double *XYZ_S,
                                 const double *normals, double *XYZ_CmS, const unsigned int BC)
{
	// Initialize DB Parameters
	char         *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int n, i, dim, FoundD, Indn;
	double       r, D, D1, D2, ABC[3], XYZ[DMAX] = {0.0}, XYZ_C[DMAX] ={0.0};

	Indn = 0;
	if (strstr(Geometry,"n-Ball")) {
		r = get_radius(Nn,(double *) XYZ_S);

		for (n = 0; n < Nn; n++) {
			if (curved_normal)
				Indn = n;

			XYZ[DMAX-1] = 0.0; // For 2D.
			for (dim = 0; dim < d; dim++)
				XYZ[dim] = XYZ_S[n+dim*Nn];

			ABC[0] = 0.0;
			ABC[1] = 0.0;
			ABC[2] = -r*r;
			for (dim = 0; dim < d; dim++) {
				ABC[0] += normals[Indn*d+dim]*normals[Indn*d+dim];
				ABC[1] += normals[Indn*d+dim]*XYZ[dim];
				ABC[2] += XYZ[dim]*XYZ[dim];
			}
			ABC[1] *= 2.0;

			// Solve quadratic equation

			// Check which solution gives the correct radius
			FoundD = 0;
			for (i = 0; i < 2; i++) {
				if (!i)
					D = (-ABC[1]+sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);
				else
					D = (-ABC[1]-sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);

				// If on opposite side of (d-1)-sphere, go to next option
				if (fabs(D) > r)
					continue;

				for (dim = 0; dim < d; dim++)
					XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;

				FoundD = 1;
				break;
			}

			if (!FoundD)
				printf("Error: Correct distance not found.\n"), EXIT_MSG;

			for (dim = 0; dim < d; dim++)
				XYZ_CmS[n+Nn*dim] = XYZ_C[dim]-XYZ[dim];
		}
	} else if (strstr(Geometry,"n-Ellipsoid")) {
		double *abc;

		abc = malloc(DMAX * sizeof *abc); // free
		get_abc_ellipse(Nn,(double *) XYZ_S,abc);

		for (n = 0; n < Nn; n++) {
			if (curved_normal)
				Indn = n;

			XYZ[DMAX-1] = 0.0; // For 2D.
			for (dim = 0; dim < d; dim++)
				XYZ[dim] = XYZ_S[n+dim*Nn];

			ABC[0] = 0.0;
			ABC[1] = 0.0;
			ABC[2] = -1.0;
			for (dim = 0; dim < d; dim++) {
				ABC[0] += normals[Indn*d+dim]*normals[Indn*d+dim]/(abc[dim]*abc[dim]);
				ABC[1] += normals[Indn*d+dim]*XYZ[dim]/(abc[dim]*abc[dim]);
				ABC[2] += XYZ[dim]*XYZ[dim]/(abc[dim]*abc[dim]);
			}
			ABC[1] *= 2.0;

			// Solve quadratic equation

			// Check which solution gives the correct radius
			r = 0.0;
			for (dim = 0; dim < d; dim++) {
				if (abc[dim] > r)
					r = abc[dim];
			}

			D1 = (-ABC[1]+sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);
			D2 = (-ABC[1]-sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);

			if (fabs(D1) < fabs(D2))
				D = D1;
			else
				D = D2;

			for (dim = 0; dim < d; dim++)
				XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;

			for (dim = 0; dim < d; dim++)
				XYZ_CmS[n+Nn*dim] = XYZ_C[dim]-XYZ[dim];
		}
		free(abc);
	} else if (strstr(Geometry,"Ringleb")) {
		if (d != 2)
			printf("Error: Unsupported.\n"), EXIT_MSG;

		unsigned int i, iMax, RinglebType = 0;
		double       Q0, KMin, KMax, q, k, a, rho, J, x, y, xS, yS, nx, ny, sign_y, alpha, beta;

		Q0   = DB.Q0;
		KMin = DB.KMin;
		KMax = DB.KMax;

		if (BC % BC_STEP_SC == BC_RIEMANN || BC % BC_STEP_SC == BC_DIRICHLET)
			RinglebType = 'f';
		else if (BC % BC_STEP_SC == BC_SLIPWALL || BC % BC_STEP_SC == BC_NEUMANN)
			RinglebType = 'w';

		sign_y = 1.0;
		if (XYZ_S[0+1*Nn] < 0.0)
			sign_y = -1.0;

		if (RinglebType == 'f') {
			q = Q0;

			for (n = 0; n < Nn; n++) {
				if (curved_normal)
					Indn = n;

				nx = normals[Indn*d+0];
				ny = normals[Indn*d+1];

				xS = XYZ_S[n     ];
				yS = XYZ_S[n+1*Nn];

				a   = sqrt(1.0-0.5*GM1*q*q);
				rho = pow(a,2.0/GM1);
				J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

				k = 0.5*(KMin+KMax);
				for (i = 0, iMax = 25; i < iMax; i++) {
					// Find point on Ringleb surface for guessed k
					x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-0.5*J;
					y = sign_y/(k*rho*q)*sqrt(1.0-pow(q/k,2.0));

					// Find the point on the normal line which is closest to the located point on the boundary
					alpha =  nx*(x-xS)+ny*(y-yS);
					beta  = -ny*(x-xS)+nx*(y-yS);

					// Find new point on the surface based on x-coordinate found above and update k
					x = xS + nx*alpha;
					k = sqrt(2.0/(2.0*rho*(x+0.5*J)+1.0/(q*q)));

//printf(" % .3e % .3e % .3e % .3e % .3e % .3e % .3e % .3e\n",x,y,xS,yS,nx,ny,alpha,beta);
					if (fabs(beta) < EPS)
						break;
				}
				if (i == iMax)
					printf("Error: Not converging (f) (% .3e).\n",beta), EXIT_MSG;

				XYZ_CmS[n+Nn*0] = x-xS;
				XYZ_CmS[n+Nn*1] = y-yS;
			}
		} else if (RinglebType == 'w') {
			double xMin,  xMax, aMin, aMax, aStep, a;

			if (XYZ_S[0+0*Nn] < 0.0) {
				k = KMax;
			} else {
				k = KMin;
			}

			xMin = xMax = XYZ_S[0];
			for (n = 1; n < Nn; n++) {
				xS = XYZ_S[n];
				if (xS < xMin) xMin = xS;
				if (xS > xMax) xMax = xS;
			}

			for (n = 0; n < Nn; n++) {
				if (curved_normal)
					Indn = n;

				nx = normals[Indn*d+0];
				ny = normals[Indn*d+1];

				xS = XYZ_S[n     ];
				yS = XYZ_S[n+1*Nn];

				// Compute normal displacement using the bisection method (ToBeModified if slow)
				aMin = (xMin-xS)/nx;
				aMax = (xMax-xS)/nx;

				aStep = fabs(0.5*(aMax-aMin));
				a = 0.5*(aMin+aMax);

				// If node is already on the boundary, set aStep = 0.0 to avoid the loop below
				x = xS;
				y = yS;
				Ringleb_boundary(&x,&y,0.0,k,'w');
				if (fabs(y-yS) < EPS)
					aStep = 0.0;
				else
					y = yS;

				while (fabs(aStep) > EPS) {
					aStep *= 0.5;

					x = xS+a*nx;
					Ringleb_boundary(&x,&y,0.0,k,'w');

					if ((y-(yS+ny*a))/ny > 0.0)
						a += aStep;
					else
						a -= aStep;
				}
				XYZ_CmS[n+Nn*0] = x-xS;
				XYZ_CmS[n+Nn*1] = y-yS;
			}
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
		}
	} else if (strstr(Geometry,"HoldenRamp")) {
		// Initialize DB Parameters
		double r = DB.rIn;

		// Standard datatypes
		unsigned int OnStraight;
		double       Xc, Yc, XcR, tan_15deg, nx, ny, x, XYZc[DMAX];

		tan_15deg = tan(15.0/180.0*PI);

		Xc = -r/tan(82.5/180.0*PI);
		Yc =  r;
		XcR = Xc + r*cos(3.0/2.0*PI+15.0/180.0*PI);

		XYZc[0] = Xc;
		XYZc[1] = Yc;
		XYZc[2] = 0.0;

		for (n = 0; n < Nn; n++) {
			if (curved_normal)
				Indn = n;

			XYZ[DMAX-1] = 0.0; // For 2D.
			for (dim = 0; dim < d; dim++)
				XYZ[dim] = XYZ_S[n+dim*Nn];

			// Find out which segment the point intersects with
			OnStraight = 0;

			// Straight line segment 1
			D = -XYZ[1]/normals[Indn*d+1];
			x =  XYZ[0]+normals[Indn*d+0]*D;
			if (x-Xc < EPS) {
				OnStraight++;
				for (dim = 0; dim < d; dim++)
					XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;
			}

			// Straight line segment 2
			nx = normals[Indn*d+0];
			ny = normals[Indn*d+1];

			// Solve for scaling of normal
			D = (XYZ[1]-tan_15deg*XYZ[0])/(-ny+tan_15deg*nx);
			x =  XYZ[0]+normals[Indn*d+0]*D;
			if (x-XcR > EPS) {
				OnStraight++;
				for (dim = 0; dim < d; dim++)
					XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;
			}

			if (OnStraight == 2)
				printf("Error: Point on both straight faces.\n"), EXIT_MSG;

			// Curve line segment
			if (!OnStraight) {
				ABC[0] = 0.0;
				ABC[1] = 0.0;
				ABC[2] = -r*r;
				for (dim = 0; dim < d; dim++) {
					ABC[0] += normals[Indn*d+dim]*normals[Indn*d+dim];
					ABC[1] += normals[Indn*d+dim]*(XYZ[dim]-XYZc[dim]);
					ABC[2] += pow(XYZ[dim]-XYZc[dim],2.0);
				}
				ABC[1] *= 2.0;

				// Solve quadratic equation

				// Check which solution gives the correct radius
				FoundD = 0;
				for (i = 0; i < 2; i++) {
					if (!i)
						D = (-ABC[1]+sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);
					else
						D = (-ABC[1]-sqrt(ABC[1]*ABC[1]-4.0*ABC[0]*ABC[2]))/(2.0*ABC[0]);

					// If on opposite side of (d-1)-sphere, go to next option
					if (fabs(D) > r)
						continue;

					FoundD = 1;
					break;
				}

				if (!FoundD)
					printf("Error: Correct distance not found.\n"), EXIT_MSG;

				for (dim = 0; dim < d; dim++)
					XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;
			}
			for (dim = 0; dim < d; dim++)
				XYZ_CmS[n+Nn*dim] = XYZ_C[dim]-XYZ[dim];
		}
	} else if (strstr(Geometry,"GaussianBump")  ||
	           strstr(Geometry,"NacaSymmetric") ||
	           strstr(Geometry,"JoukowskiSymmetric") ||
	           strstr(Geometry,"ExpansionCorner") ||
	           strstr(Geometry,"EllipsoidalBump")) {
		double       nx, ny, xS, yS, DStep, x, y, h, xp, yp, ypE;

		unsigned int i, j, count, countMax;
		double       Dmult, y_ratio, yE, Dup;

		surface_tdef f_surface;

		select_functions_surface(&f_surface);

		if (d == 3)
			printf("Add support.\n"), EXIT_MSG;

		// Obtain approximate measure of mesh size
		// Note: Distance between straight and curved surface ~= O(h^2)
		h = 0.0;
		for (i = 0;   i < Nn; i++) {
		for (j = i+1; j < Nn; j++) {
			D = 0.0;
			for (dim = 0; dim < d; dim++)
				D += pow(XYZ_S[i+dim*Nn]-XYZ_S[j+dim*Nn],2.0);
			D = sqrt(D);
			if (D > h)
				h = D;
		}}

		for (n = 0; n < Nn; n++) {
			if (curved_normal)
				Indn = n;

			nx = normals[Indn*d+0];
			ny = normals[Indn*d+1];

			XYZ[DMAX-1] = 0.0; // For 2D.
			for (dim = 0; dim < d; dim++)
				XYZ[dim] = XYZ_S[n+dim*Nn];

			// Use the Bisection Method to solve for D.
			// Note: Newton's method is slow because of nearly zero slope in the gaussian.
			xS = XYZ_S[n     ];
			yS = XYZ_S[n+1*Nn];

			Dmult = -1.0;
			// Check if the normal points towards or away from the surface
			x = xS;
			y = yS;

			if (d == 2)
				D = yS-f_surface(x,y,d);
			else
				printf("Add support\n"), EXIT_MSG;

			if (D/ny < 0.0)
				Dmult = 1.0;

			// Find Distance from node to surface
			if (h < 1.0)
				DStep = 5e-2*Dmult*h*h;
			else
				DStep = 1e-1*Dmult;

			// Check if the node is on the surface
			if (fabs(D) < EPS)
				DStep = 0.0;

			D = 0.0;
			if (fabs(DStep) > EPS) {
				// Step towards the surface with step DStep until close enough to ensure that it can be reached using
				// the bisection method below.
				countMax = 100;
				for (count = 0; count < countMax; count++) {
					Dup = D+DStep;

					xp  = xS+D*nx;
					yp  = yS+D*ny;
					ypE = f_surface(xp,yp,d);

					x  = xS+Dup*nx;
					y  = yS+Dup*ny;
					yE = f_surface(x,y,d);

					y_ratio = (yp-ypE)/(y-yE);

					if ((fabs(yp-ypE) < fabs(y-yE)) || y_ratio < 0.0)
						break;

					D = Dup;
				}
				if (count == countMax)
					printf("Error: Potential problem (% .3e % .3e).\n",yp-ypE,y-yE), EXIT_MSG;
			}

			while (fabs(DStep) > 1e-1*EPS) {
				Dup = D+DStep;

				xp  = xS+D*nx;
				yp  = yS+D*ny;
				ypE = f_surface(xp,yp,d);

				x  = xS+Dup*nx;
				y  = yS+Dup*ny;
				yE = f_surface(x,y,d);

				y_ratio = (yp-ypE)/(y-yE);

				if ((fabs(yp-ypE) > fabs(y-yE)) && y_ratio > 0.0)
					D = Dup;

				DStep *= 0.5;
			}

			// Check that the surface was reached
			x = xS+D*nx;
			y = yS+D*ny;

			if (d == 2) {
				yE = f_surface(x,y,d);
//				if (fabs(y-yE) > 1e2*EPS) {
				if (fabs(y-yE) > 2e1*SQRT_EPS) {
					printf("% .3e % .3e\n",xS,yS);
					printf("Error: Did not reach the surface (% .3e % .3e % .3e % .3e % .3e).\n",
					       D,x,y,yE,y-yE), EXIT_MSG;
				}
			} else {
				printf("Add support.\n"), EXIT_MSG;
			}

			for (dim = 0; dim < d; dim++) {
				XYZ_C[dim] = XYZ[dim]+normals[Indn*d+dim]*D;
				XYZ_CmS[n+Nn*dim] = XYZ_C[dim]-XYZ[dim];
			}
		}
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

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

			if (Z-rXYZ > 0.0)
				p = acos(1.0);
			else
				p = acos(Z/rXYZ);
		}

		PComps[0*Nn+n] = t;
		PComps[1*Nn+n] = p;
	}

	if (d == 3) {
		// If one (and only one) of the vertices lies on a pole of the sphere, re-calculate theta for the vertex based
		// on the theta values of the other vertices.
		for (n = 0; n < Nn; n++) {
			if (VeXYZ[n][2]-r > 0.0)
				p = acos(1.0);
			else
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
		} else {
			printf("Error: Unsupported.\n"), EXIT_MSG;
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

static void compute_pc_ellipsoid(struct S_pc *data) {
	/*
	 *	Purpose:
	 *		Compute (t)heta and (p)hi parametrization components associated with the input vertices.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn;
	double       a, b, t, p, *PComps, **VeXYZ, X, Y;

	if (d == 3)
		printf("Add support.\n"), EXIT_MSG;

	if (data->VeSurface == 0) {
		a = DB.aIn;
		b = DB.bIn;
	} else if (data->VeSurface == 1) {
		a = DB.aOut;
		b = DB.bOut;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	Nn     = data->Nn;
	VeXYZ  = data->VeXYZ;
	PComps = data->PComps;

	for (n = 0; n < Nn; n++) {
		X = VeXYZ[n][0];
		Y = VeXYZ[n][1];

		t = atan2(Y*a,X*b);

		if (d == 2) {
			p = PI/2.0;
		} else {
			printf("Add support.\n"), EXIT_MSG;
		}

		PComps[0*Nn+n] = t;
		PComps[1*Nn+n] = p;
	}
}

static void compute_XYZ_ellipsoid(struct S_XYZ *data)
{
	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, Nn;
	double       *theta, *phi, *XYZ, *X, *Y, *Z, t, p, a, b, c;

	if (d == 3)
		printf("Add support.\n"), EXIT_MSG;

	if (data->VeSurface == 0) {
		a = DB.aIn;
		b = DB.bIn;
		c = DB.cIn;
	} else if (data->VeSurface == 1) {
		a = DB.aOut;
		b = DB.bOut;
		c = DB.cOut;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

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

		X[n] = a*cos(t)*sin(p);
		Y[n] = b*sin(t)*sin(p);
		if (d == 3)
			Z[n] = c*cos(p);
	}
}

static void select_functions_Curved(compute_pc_tdef *compute_pc, compute_XYZ_tdef *compute_XYZ)
{
	// Initialize DB Parameters
	char *Geometry = DB.Geometry;

	if (strstr(Geometry,"n-Ball")) {
		*compute_pc  = compute_pc_dsphere;
		*compute_XYZ = compute_XYZ_dsphere;
	} else if (strstr(Geometry,"n-Ellipsoid")) {
		*compute_pc  = compute_pc_ellipsoid;
		*compute_XYZ = compute_XYZ_ellipsoid;
		if (DB.Parametrization != NORMAL &&
		    DB.Parametrization != RADIAL_PROJECTION)
				printf("Error: Unsupported parametrization.\n"), EXIT_MSG;
	} else if (strstr(Geometry,"Ringleb")       ||
			   strstr(Geometry,"HoldenRamp")    ||
			   strstr(Geometry,"GaussianBump")  ||
			   strstr(Geometry,"NacaSymmetric") ||
			   strstr(Geometry,"JoukowskiSymmetric") ||
			   strstr(Geometry,"ExpansionCorner") ||
			   strstr(Geometry,"EllipsoidalBump")) {
		if (DB.Parametrization != NORMAL)
			printf("Add support if not using NORMAL parametrization.\n"), EXIT_MSG;
		*compute_pc  = NULL;
		*compute_XYZ = NULL;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

static double *compute_BlendV(struct S_Blend *data, const unsigned int order)
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
	double       *I_vGs_vGc, *I_bGs_vGc, *BlendV, BlendNum, BlendDen;

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

	if (Blending == SCOTT && type == TRI) {
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
			BlendV[n] = pow(1.0-I_vGs_vGc[n*Nve+b],(double) order);
		}
	} else if (EclassV == C_SI) { // Blending == SZABO_BABUSKA
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
	} else if (EclassV == C_TP) { // Blending == GORDON_HALL
		for (n = 0; n < NvnG; n++) {
			BlendV[n] = 0.0;
			for (ve = 0; ve < Nbve[b]; ve++)
				BlendV[n] += I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]];
		}
	} else {
		printf("%d %d\n",Blending,EclassV);
		printf("Error: Unsupported.\n"), EXIT_MSG;
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
	unsigned int b, n, ve, dim, NvnG, NbnG, Nve, *Nbve, *VeBcon, *VeInfo, NbveMax, BC;
	double       **VeXYZ, *XYZ, *XYZ_vV, *PComps_V, *PComps_B, *XYZ_CmS, *XYZ_update, *XYZ_S, *n_S, d_S, nNorm, r,
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
	Nbve    = data->Nbve;
	VeInfo  = data->VeInfo;
	VeBcon  = data->VeBcon;
	NbveMax = data->NbveMax;
	BC      = data->BC;

	XYZ       = data->XYZ;
	XYZ_vV    = data->XYZ_vV;
	I_vGs_vGc = data->I_vGs_vGc;
	I_vGc_bGc = data->I_vGc_bGc;
	I_bGc_vGc = data->I_bGc_vGc;

	data_pc  = malloc(sizeof *data_pc);  // free
	data_XYZ = malloc(sizeof *data_XYZ); // free

	compute_pc  = data->compute_pc;
	compute_XYZ = data->compute_XYZ;

	if (Parametrization == ARC_LENGTH) {
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

		// Subtract XYZ_S;
		mm_d(CBCM,CBT,CBNT,NbnG,d,NvnG,-1.0,1.0,I_vGc_bGc,XYZ,XYZ_CmS);
	} else if (Parametrization == RADIAL_PROJECTION) {
		// Compute XYZ_C
		XYZ_S = malloc(NbnG*d * sizeof *XYZ_S); // free

		mm_d(CBCM,CBT,CBNT,NbnG,d,NvnG,1.0,0.0,I_vGc_bGc,XYZ,XYZ_S);

		VeXYZ = malloc(NbnG * sizeof *VeXYZ); // free
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
	} else if (Parametrization == NORMAL || Parametrization == ORDER_H) {
// Change comment here as this may no longer be XYZ_S (ToBeDeleted)
		// Compute XYZ_S
		XYZ_S = malloc(NbnG*d * sizeof *XYZ_S); // free
		mm_d(CBCM,CBT,CBNT,NbnG,d,NvnG,1.0,0.0,I_vGc_bGc,XYZ,XYZ_S);

		// Find coordinates of 3 vertices used to compute the normal to the boundary
		VeXYZ = malloc(DMAX * sizeof *VeXYZ); // free
		if (d == 3) {
			for (ve = 0; ve < DMAX; ve++) {
				VeXYZ[ve] = calloc(DMAX , sizeof **VeXYZ); // free
				for (dim = 0; dim < d; dim++)
					VeXYZ[ve][dim] = XYZ_vV[Nve*dim+VeBcon[b*NbveMax+ve]];
			}
		} else if (d == 2) {
			for (ve = 0; ve < d; ve++) {
				VeXYZ[ve] = calloc(DMAX , sizeof **VeXYZ); // free
				for (dim = 0; dim < d; dim++)
					VeXYZ[ve][dim] = XYZ_vV[Nve*dim+VeBcon[b*NbveMax+ve]];
			}
			ve = d;
			VeXYZ[ve] = calloc(DMAX , sizeof **VeXYZ); // free
			for (dim = 0; dim < d; dim++)
				VeXYZ[ve][dim] = XYZ_vV[Nve*dim+VeBcon[b*NbveMax+0]];
			// Set z component to something non-zero
			VeXYZ[ve][d] = 1.0;
		}

		n_S = malloc(DMAX * sizeof *n_S); // free
		// Compute normal vector to the straight FACE
		compute_plane(VeXYZ[0],VeXYZ[1],VeXYZ[2],n_S,&d_S);

		// Normalize
		nNorm = array_norm_d(d,n_S,"L2");
		for (dim = 0; dim < d; dim++)
			n_S[dim] = n_S[dim]/nNorm;
		d_S /= nNorm;

		// Ensure that normal points outwards
		r = -d_S;
		for (dim = 0; dim < d; dim++)
			r += n_S[dim]*VeXYZ[0][dim];

		if (r > 0.0) {
			for (dim = 0; dim < d; dim++)
				n_S[dim] *= -1.0;
			d_S *= -1.0;
		}

		// Compute normal distance from polynomial FACE to curved geometry
		if (Parametrization == NORMAL) {
			compute_normal_displacement(NbnG,0,XYZ_S,n_S,XYZ_CmS,BC);
		} else if (Parametrization == ORDER_H) {
			// Compute XYZ_CmS using XYZ_S + perturbation
			unsigned int Blending = DB.Blending;

			double h, r, *XYZ_Sp, *VeXYZb;

			// Compute h (surface length)
			h = 0.0;
			if (d == 2) {
				for (dim = 0; dim < d; dim++)
					h += pow(VeXYZ[0][dim]-VeXYZ[1][dim],2.0);
				h = sqrt(h);
			} else {
				printf("Error: Unsupported.\n"), EXIT_MSG;
			}

			// Compute perturbed XYZ_S nodes

			// Find barycentric coordinates relating straight and curved BOUNDARY geometry nodes
			I_vGs_bGc = mm_Alloc_d(CBRM,CBNT,CBNT,NbnG,Nve,NvnG,1.0,I_vGc_bGc,I_vGs_vGc); // free

			I_bGs_bGc = malloc(NbnG*Nbve[b] * sizeof *I_bGs_bGc); // free
			for (n = 0; n < NbnG; n++) {
			for (ve = 0; ve < Nbve[b]; ve++) {
				I_bGs_bGc[n*Nbve[b]+ve] = I_vGs_bGc[n*Nve+VeBcon[b*NbveMax+ve]];
			}}
			free(I_vGs_bGc);

			// Perturb the coordinates of I_bGs_bGc
			if (Blending != SZABO_BABUSKA || d != 2)
				printf("Error: Unsupported.\n"), EXIT_MSG;

			for (n = 0; n < NbnG; n++) {
				r = I_bGs_bGc[n*Nbve[b]+1]-I_bGs_bGc[n*Nbve[b]+0];
				if (fabs(fabs(r)-1.0) > EPS) {
					// A.5 refers to equation A.5 in Zwanenburg(2017)-ToBeModified
					r *= (1+0.125*(1-r*r)*pow(h,2.0)); // ~Radial (Truncated after 1st term, A.5 (i = 2, C_1 = 1/8))
//					r *= (1-(1-r*r)*pow(h,2.0));       // OK      (A.5 (i = 2, C_1 = -1))
//					r *= (1-(1-r*r)*pow(h,1.0));       // NOT OK  (A.5 (i = 1))

					I_bGs_bGc[n*Nbve[b]+0] = 0.5*(1.0-r);
					I_bGs_bGc[n*Nbve[b]+1] = 0.5*(1.0+r);
				}
			}

			// Compute perturbed BOUNDARY nodes (XYZ_Sp)
			VeXYZb = malloc(Nbve[b]*d * sizeof *VeXYZb); // free
			for (ve = 0; ve < Nbve[b]; ve++) {
			for (dim = 0; dim < d; dim++) {
				VeXYZb[dim*Nbve[b]+ve] = VeXYZ[ve][dim];
			}}

			XYZ_Sp = malloc(NbnG*d * sizeof *XYZ_Sp); // free

			mm_d(CBCM,CBT,CBNT,NbnG,d,Nbve[b],1.0,0.0,I_bGs_bGc,VeXYZb,XYZ_Sp);
			free(I_bGs_bGc);
			free(VeXYZb);

			// Compute XYZ_CmS using XYZ_Sp
			compute_normal_displacement(NbnG,0,XYZ_Sp,n_S,XYZ_CmS,BC);

			// Correct XYZ_CmS so that straight contribution is from XYZ_S
			for (n = 0; n < NbnG*d; n++)
				XYZ_CmS[n] += XYZ_Sp[n]-XYZ_S[n];

			free(XYZ_Sp);
		}
		array_free2_d(DMAX,VeXYZ);

		// Set displacement of previously existing vertices to 0.0 (Needed for Ringleb)
		double XYZdiff;
		for (ve = 0; ve < Nve; ve++) {
			for (n = 0; n < NbnG; n++) {
				XYZdiff = 0.0;
				for (dim = 0; dim < d; dim++)
					XYZdiff += pow(XYZ_S[dim*NbnG+n]-XYZ_vV[dim*Nve+VeBcon[b*NbveMax+ve]],2.0);
				XYZdiff = sqrt(XYZdiff);

				if (XYZdiff < EPS) {
					for (dim = 0; dim < d; dim++)
						XYZ_CmS[dim*NbnG+n] = 0.0;
					break;
				}
			}
		}

		free(n_S);
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

static void blend_boundary(struct S_VOLUME *VOLUME, const unsigned int BType, const unsigned int vertex_blending)
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
	unsigned int d           = DB.d,
	             Blending_HO = DB.Blending_HO;

	// Standard datatypes
	unsigned int b, dim, n, Vb, P, PV, PMin, IndSurf, Vtype,
	             Nve, Nb, *Nbve, *VeBcon, NbveMax, NbrefMax, NvnG, BC, *VeInfo;
	double       *XYZ, *XYZ_update, *XYZ_vV, *BlendV, *I_vGs_vGc;

	struct S_ELEMENT *ELEMENT, *ELEMENT_B;
	struct S_pc      *data_pc;
	struct S_XYZ     *data_XYZ;
	struct S_Blend   *data_blend;

	// silence
	Nb = NbrefMax = NbveMax = IndSurf = BC = 0;
	VeBcon = Nbve = NULL;
	ELEMENT_B = NULL;

	data_pc    = malloc(sizeof *data_pc);    // free
	data_XYZ   = malloc(sizeof *data_XYZ);   // free
	data_blend = malloc(sizeof *data_blend); // free

	select_functions_Curved(&data_blend->compute_pc,&data_blend->compute_XYZ);

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	if (vertex_blending) {
		PV = 2;

		NvnG = ELEMENT->NveP2;
		XYZ  = VOLUME->XYZ_vVP2;
	} else {
		PV   = VOLUME->P;
		NvnG = VOLUME->NvnG;
		XYZ  = VOLUME->XYZ;
	}

	Vtype  = VOLUME->type;
	XYZ_vV = VOLUME->XYZ_vV;
	VeInfo = VOLUME->VeInfo;

	data_blend->NvnG   = NvnG;
	data_blend->XYZ    = XYZ;
	data_blend->XYZ_vV = XYZ_vV;
	data_blend->VeInfo = VeInfo;

	data_blend->type = Vtype;

	Nve       = ELEMENT->Nve;
	if (vertex_blending)
		I_vGs_vGc = ELEMENT->I_vGs_vG2[1][PV][0];
	else
		I_vGs_vGc = ELEMENT->I_vGs_vGc[1][PV][0];

	data_blend->EclassV   = ELEMENT->Eclass;
	data_blend->Nve       = Nve;
	data_blend->I_vGs_vGc = I_vGs_vGc;

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

	if (!Blending_HO)
		PMin = PV;
	else
		PMin = 2;

	for (P = PMin; P <= PV; P++) { // Note: For vertex_blending = 1, PMin == PV == 2
	for (b = 0; b < Nb; b++) {
		if (BType == 'e')
			BC = VOLUME->BC[1][b];
		else if (BType == 'f')
			BC = VOLUME->BC[0][b];

		if (BC / BC_STEP_SC != 2)
			continue;

		Vb = b*NbrefMax;
		if (BType == 'e') {
			ELEMENT_B = get_ELEMENT_type(LINE);
			if (vertex_blending) {
				data_blend->I_vGc_bGc = ELEMENT->I_vG2_eG2[PV][P][Vb];
				data_blend->I_bGc_vGc = ELEMENT->I_eG2_vG2[P][PV][Vb];
				data_blend->I_bGs_vGc = ELEMENT->I_eGs_vG2[1][PV][Vb];
			} else {
				data_blend->I_vGc_bGc = ELEMENT->I_vGc_eGc[PV][P][Vb];
				data_blend->I_bGc_vGc = ELEMENT->I_eGc_vGc[P][PV][Vb];
				data_blend->I_bGs_vGc = ELEMENT->I_eGs_vGc[1][PV][Vb];
			}
		} else if (BType == 'f') {
			ELEMENT_B = get_ELEMENT_F_type(ELEMENT->type,b);
			if (vertex_blending) {
				data_blend->I_vGc_bGc = ELEMENT->I_vG2_fG2[PV][P][Vb];
				data_blend->I_bGc_vGc = ELEMENT->I_fG2_vG2[P][PV][Vb];
				data_blend->I_bGs_vGc = ELEMENT->I_fGs_vG2[1][PV][Vb];
			} else {
				data_blend->I_vGc_bGc = ELEMENT->I_vGc_fGc[PV][P][Vb];
				data_blend->I_bGc_vGc = ELEMENT->I_fGc_vGc[P][PV][Vb];
				data_blend->I_bGs_vGc = ELEMENT->I_fGs_vGc[1][PV][Vb];
			}
		}

		if (vertex_blending)
			data_blend->NbnG = ELEMENT_B->NvnG2[P];
		else
			data_blend->NbnG = ELEMENT_B->NvnGc[P];

		data_blend->b    = b;
		data_blend->BC   = BC;

		if (!Blending_HO) {
			BlendV = compute_BlendV(data_blend,2);   // free
		} else {
			BlendV = compute_BlendV(data_blend,P);   // free
		}
		XYZ_update = compute_XYZ_update(data_blend); // free

		// Blend BOUNDARY perturbation to VOLUME geometry nodes
		for (n = 0; n < NvnG; n++) {
			// Don't move internal VOLUME nodes when moving vertices to the exact geometry.
			if ((BlendV[n] < EPS) || (fabs(BlendV[n]-1.0) > EPS && vertex_blending))
				continue;

			for (dim = 0; dim < d; dim++)
				XYZ[n+NvnG*dim] += BlendV[n]*XYZ_update[n+NvnG*dim];
		}

		free(XYZ_update);
		free(BlendV);
	}}

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
		blend_boundary(VOLUME,'e',0);

	// Treat curved FACEs
	if (Vcurved == 1)
		blend_boundary(VOLUME,'f',0);
}

void setup_Curved_vertices(struct S_VOLUME *VOLUME)
{
	/*
	 *	Purpose:
	 *		Ensure that newly created vertices are placed on the boundary, in well behaved positions.
	 *
	 *	Comments:
	 *		See comments in 'vertices_to_exact_geom' for why blending is used here. For example, for the Ringleb case,
	 *		computing the new vertex location assuming correct x-coordinate and computing the y-coordinate results in an
	 *		extremely stretched mesh near the y = 0 axis.
	 */

	// Initialize DB Parameters
	unsigned int d    = DB.d;

	// Standard datatypes
	unsigned int Nve, NveP2, Vcurved;
	double       *XYZ_vV, *XYZ_vVP2, *I_vGs_vGsP2;

	struct S_ELEMENT *ELEMENT;

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	Nve   = ELEMENT->Nve;
	NveP2 = ELEMENT->NveP2;

	I_vGs_vGsP2 = ELEMENT->I_vGs_vGs[1][2][0];

	if (VOLUME->XYZ_vVP2)
		free(VOLUME->XYZ_vVP2);
	VOLUME->XYZ_vVP2 = malloc(NveP2*d * sizeof *(VOLUME->XYZ_vVP2)); // keep

	// Compute straight geometry nodes of P2 element
	XYZ_vV   = VOLUME->XYZ_vV;
	XYZ_vVP2 = VOLUME->XYZ_vVP2;
	mm_CTN_d(NveP2,d,Nve,I_vGs_vGsP2,XYZ_vV,XYZ_vVP2);

	Vcurved = VOLUME->curved;

	// Treat curved EDGEs not on curved FACEs (3D only)
	if (d == DMAX)
		blend_boundary(VOLUME,'e',1);

	// Treat curved FACEs
	if (Vcurved == 1)
		blend_boundary(VOLUME,'f',1);
}
