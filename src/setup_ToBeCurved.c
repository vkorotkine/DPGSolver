// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "setup_ToBeCurved.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "array_norm.h"
#include "cubature.h"
#include "element_functions.h"
#include "matrix_functions.h"

#include "array_print.h" // ToBeDeleted

/*
 *	Purpose:
 *		Set up "ToBeCurved" meshes.
 *
 *	Comments:
 *		Requires that straight VOLUME nodes have already been stored.
 *		May want to change parameters in PeriodicVortex case to take the scaling into account (ToBeDeleted). The orders
 *		obtained seem to be correct however.
 *
 *	Notation:
 *		XYZ   : High-order VOLUME (XYZ) coordinates after curving.
 *
 *	References:
 *		Gordon(1973)-Construction_of_Curvilinear_Coordinate_Systems_and_Applications_to_Mesh_Generation
 *		Gordon(1973)-Transfinite_Element_Methods-Blending-Function_Interpolation_over_Arbitraty_Curved_Element_Domains
 *		Rosca(2011)-Uniform_Spherical_Grids_via_Equal_Area_Projection_from_the_Cube_to_the_Sphere
 */

static void         ToBeCurved_elliptic_pipe    (unsigned int const Nn, double const *const XYZ_S, double *const XYZ);
static void         ToBeCurved_parabolic_pipe    (unsigned int const Nn, double const *const XYZ_S, double *const XYZ);
static void         ToBeCurved_sinusoidal_pipe    (unsigned int const Nn, double const *const XYZ_S, double *const XYZ);
static void         ToBeCurved_cube_to_sphere   (unsigned int Nn, double *XYZ_S, double *XYZ);
static void         ToBeCurved_square_to_circle (unsigned int Nn, double *XYZ_S, double *XYZ);
static double         *cube_to_sphere           (double XY[2], unsigned int OrderOut[3], int SignOut, double beta);
static void         ToBeCurved_TP               (unsigned int Nn, double *XYZ_S, double *XYZ);
static unsigned int   is_in_blend_region        (double XYZn[3]);
static void           get_blend_bounds          (const double Xn, const double Zn, unsigned int *xLoc, unsigned int *zLoc,
                                                 double *xBounds, double *zBounds);
static double         get_arc_length            (const double XL, const double XR, const double Z,
                                                 const unsigned int DOrder[2]);
static double         *eval_TP_function         (const unsigned int Nn, const double *XZ, const unsigned int DOrder[2],
                                                 const unsigned int Single, double **abcP);

static void ToBeCurved_elliptic_pipe (unsigned int const Nn, double const *const XYZ_S, double *const XYZ)
{
	if (DB.d != 2)
		EXIT_UNSUPPORTED;

	double *const X = &XYZ[Nn*0],
	       *const Y = &XYZ[Nn*1];

	double const *const X_S = &XYZ_S[Nn*0],
	             *const Y_S = &XYZ_S[Nn*1];

	double const a = 2.0, b = 4.0, a_1 = 1.0, a_2 = 2.0;


	for (size_t n = 0; n < Nn; n++) {
	     Y[n] = (b/(2*a))*((a_2-a_1)*Y_S[n]+(a_2+a_1))*sin((PI/2)*(X_S[n]+1));
             X[n] = -0.5*((a_2-a_1)*Y_S[n]+(a_2+a_1))*cos((PI/2)*(X_S[n]+1));
	}
}

static void ToBeCurved_parabolic_pipe (unsigned int const Nn, double const *const XYZ_S, double *const XYZ)
{
        if (DB.d != 2)
                EXIT_UNSUPPORTED;

        double *const X = &XYZ[Nn*0],
               *const Y = &XYZ[Nn*1];

        double const *const X_S = &XYZ_S[Nn*0],
                     *const Y_S = &XYZ_S[Nn*1];

        double const  b = 0.25, a_1 = 0.25, a_2 = 2.25;

        for (size_t n = 0; n < Nn; n++) {

			 Y[n] = (a_2+a_1)/2+(a_2-a_1)*Y_S[n]/2-a_1*pow((X_S[n]+1),2)/4;
             X[n] = sqrt(a_1/b)*(X_S[n]+1)/2;
        }
}

static void ToBeCurved_sinusoidal_pipe (unsigned int const Nn, double const *const XYZ_S, double *const XYZ)
{
        if (DB.d != 2)
                EXIT_UNSUPPORTED;

        double *const X = &XYZ[Nn*0],
               *const Y = &XYZ[Nn*1];

        double const *const X_S = &XYZ_S[Nn*0],
                     *const Y_S = &XYZ_S[Nn*1];

        double const a = 2.0,  b = PI/2, c_1 = 0, c_2 = 2.0;

        for (size_t n = 0; n < Nn; n++) {
               Y[n] = (c_2+c_1)/2+(c_2-c_1)*Y_S[n]/2+a*cos(X_S[n]*acos(-c_1/a));
               X[n] = (1/b)*acos(-c_1/a)*X_S[n];
        }
}


static void ToBeCurved_sphere_to_ellipsoid(const unsigned int Nn, double *XYZ)
{
	/*
	 *	Comments:
	 *		x = a*cos(t)*sin(p)
	 *		y = b*sin(t)*sin(p)
	 *		x = c*cos(p)
	 */

	// Initialize DB Parameters
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
	double       *X, *Y, *Z, r2, r, t, p, ratio, a, b, c;

	X = &XYZ[Nn*0];
	Y = &XYZ[Nn*1];
	Z = &XYZ[Nn*(d-1)];

//	printf("Make modifications so that this works with high aspect ratio ellipsoids.\n"), EXIT_MSG;

	for (n = 0; n < Nn; n++) {
		r2 = X[n]*X[n]+Y[n]*Y[n];
		if (d == 3)
			r2 += Z[n]*Z[n];

		r = sqrt(r2);
		ratio = (r-rIn)/(rOut-rIn);

		a = aIn+ratio*(aOut-aIn);
		b = bIn+ratio*(bOut-bIn);
		c = cIn+ratio*(cOut-cIn);

//		t = atan2(Y[n]/b,X[n]/a); // Decreased mesh regularity for this option causes loss of optimal orders.
		t = atan2(Y[n],X[n]);

		if (d == 2) {
			p = PI/2.0;
		} else {
			if (Z[n]-c > 1.0)
				p = acos(1.0);
			else
				p = acos(Z[n]/c);
		}

		X[n] = a*cos(t)*sin(p);
		Y[n] = b*sin(t)*sin(p);
		if (d == 3)
			Z[n] = c*cos(p);
	}
}

static void correct_ToBeCurved(struct S_VOLUME *VOLUME)
{
	/*
	 *	Purpose:
	 *		Correct position of internal VOLUME geometry nodes using blending.
	 *
	 *	Comments:
	 *		This only needs to loop over element faces (even in 3D) as all curved surface nodes are given by the
	 *		ToBeCurved projection.
	 *		Only SB blending is supported for SI elements here as a result of all blending functions being equivalent.
	 */

	if (DB.Blending_HO)
		printf("Error: Unsupported.\n"), EXIT_MSG;

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, dim, ve, P, PV, f, Vf, NvnG, NfnG, Eclass, Nf, Nve, BC, *Nfve, *VeFcon, InternalCurved;
	double       *XYZ, *XYZ_S, *XYZ_C, *XYZ_CmS, *XYZ_update, *I_vGs_vGc, *I_vGc_fGc, *I_fGc_vGc, *I_fGs_vGc,
	             BlendNum, BlendDen, *BlendV;

	struct S_ELEMENT *ELEMENT, *ELEMENT_F;

	InternalCurved = 1;

	PV     = VOLUME->P;
	NvnG   = VOLUME->NvnG;
	XYZ_C  = VOLUME->XYZ;
	Eclass = VOLUME->Eclass;

	if (Eclass == C_SI && DB.Blending != SZABO_BABUSKA)
		printf("Error: Please select SB blending.\n"), EXIT_MSG;

	ELEMENT = get_ELEMENT_type(VOLUME->type);

	Nf     = ELEMENT->Nf;
	Nve    = ELEMENT->Nve;
	Nfve   = ELEMENT->Nfve;
	VeFcon = ELEMENT->VeFcon;

	I_vGs_vGc = ELEMENT->I_vGs_vGc[1][PV][0];

	// Initialize XYZ using projected vertices
	XYZ = malloc(NvnG*d * sizeof *XYZ); // VOLUME->XYZ freed then XYZ assigned
	mm_d(CBCM,CBT,CBNT,NvnG,d,Nve,1.0,0.0,ELEMENT->I_vGs_vGc[1][PV][0],VOLUME->XYZ_vVc,XYZ);

	P = PV;
	for (f = 0; f < Nf; f++) {
		if (!InternalCurved) {
			BC = VOLUME->BC[0][f];
			if (BC / BC_STEP_SC != 2)
				continue;
		}

		Vf = f*NFREFMAX;

		ELEMENT_F = get_ELEMENT_F_type(ELEMENT->type,f);

		NfnG      = ELEMENT_F->NvnGc[P];
		I_vGc_fGc = ELEMENT->I_vGc_fGc[PV][P][Vf];
		I_fGc_vGc = ELEMENT->I_fGc_vGc[P][PV][Vf];

		// Compute XYZ_update
		XYZ_CmS = malloc(NfnG*d * sizeof *XYZ_S); // free
		mm_d(CBCM,CBT,CBNT,NfnG,d,NvnG, 1.0,0.0,I_vGc_fGc,XYZ_C,XYZ_CmS);
		mm_d(CBCM,CBT,CBNT,NfnG,d,NvnG,-1.0,1.0,I_vGc_fGc,XYZ,  XYZ_CmS);

		XYZ_update = mm_Alloc_d(CBCM,CBT,CBNT,NvnG,d,NfnG,1.0,I_fGc_vGc,XYZ_CmS); // free
		free(XYZ_CmS);

		// Compute blending function
		I_fGs_vGc = ELEMENT->I_fGs_vGc[1][PV][Vf];

		BlendV = malloc(NvnG * sizeof *BlendV); // free
		if (Eclass == C_SI) { // SB Blending
			for (n = 0; n < NvnG; n++) {
				BlendNum = 1.0;
				BlendDen = 1.0;
				for (ve = 0; ve < Nfve[f]; ve++) {
					BlendNum *= I_vGs_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
					BlendDen *= I_fGs_vGc[n*Nfve[f]+ve];
				}
				if (BlendNum < EPS)
					BlendV[n] = 0.0;
				else
					BlendV[n] = BlendNum/BlendDen;
			}
		} else if (Eclass == C_TP) { // GH Blending
			for (n = 0; n < NvnG; n++) {
				BlendV[n] = 0.0;
				for (ve = 0; ve < Nfve[f]; ve++)
					BlendV[n] += I_vGs_vGc[n*Nve+VeFcon[f*NFVEMAX+ve]];
			}
		} else {
			printf("Add support.\n"), EXIT_MSG;
		}

		// Blend surface geometry perturbation to the VOLUME
		for (n = 0; n < NvnG; n++) {
			if (BlendV[n] < EPS)
				continue;

			for (dim = 0; dim < d; dim++)
				XYZ[n+NvnG*dim] += BlendV[n]*XYZ_update[n+NvnG*dim];
		}
		free(XYZ_update);
		free(BlendV);
	}

	free(VOLUME->XYZ);
	VOLUME->XYZ = XYZ;
}

void setup_ToBeCurved(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
	             *Geometry = DB.Geometry;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int returnStraight = 0; // return the straight mesh
	unsigned int i, nG, dim, NvnG, correctTBC;
	double *XYZ, *XYZ_S;

	struct S_ELEMENT *ELEMENT;

	correctTBC = 1; // Correct internal nodes with blending
	for (nG = 0; nG < 2; nG++) {
		if (nG == 0) {
			NvnG = VOLUME->NvnG;
			XYZ_S = VOLUME->XYZ_S;

			XYZ = malloc (NvnG*d * sizeof *XYZ); // keep
			VOLUME->XYZ = XYZ;
		} else if (nG == 1) {
			ELEMENT = get_ELEMENT_type(VOLUME->type);
			NvnG = ELEMENT->Nve;
			XYZ_S = VOLUME->XYZ_vV;

			XYZ = malloc (NvnG*d * sizeof *XYZ); // keep
			VOLUME->XYZ_vVc = XYZ;
		}

		if (returnStraight) {
			for (i = 0; i < NvnG*d; i++)
				XYZ[i] = XYZ_S[i];
		} else {
			if (strstr(Geometry,"dm1-Spherical_Section")) {
					ToBeCurved_cube_to_sphere(NvnG,XYZ_S,XYZ);
			} else if (strstr(Geometry,"Ellipsoidal_Section")) {
					ToBeCurved_cube_to_sphere(NvnG,XYZ_S,XYZ);
					ToBeCurved_sphere_to_ellipsoid(NvnG,XYZ);
			} else if (strstr(TestCase,"GaussianBump") ||
					   strstr(TestCase,"PolynomialBump")) {
					ToBeCurved_TP(NvnG,XYZ_S,XYZ);
			} else if (strstr(TestCase,"PeriodicVortex")) {
				double n = 2.0, A = 0.1, L0 = 2.0, dxyz = 1.0, scale, *X0, *Y0, *Z0;

				// silence
				X0 = Y0 = Z0 = NULL;

				scale = 2.0;
				DB.PeriodL = scale*2.0;

				// d > 1 for this case
				for (dim = 0; dim < d; dim++) {
					X0 = &XYZ_S[0*NvnG];
					Y0 = &XYZ_S[1*NvnG];
					if (dim == 2)
						Z0 = &XYZ_S[2*NvnG];
				}

				if (d == 2) {
					for (i = 0; i < NvnG; i++) {
						XYZ[       i] = scale*(X0[i] + A*dxyz*sin(n*PI/L0*Y0[i]));
						XYZ[1*NvnG+i] = scale*(Y0[i] + A*dxyz*sin(n*PI/L0*X0[i]));
					}
				} else if (d == 3) {
					for (i = 0; i < NvnG; i++) {
//						XYZ[       i] = scale*X0[i];// + A*dxyz*sin(n*PI/L0*Y0[i])*sin(n*PI/L0*Z0[i]);
//						XYZ[1*NvnG+i] = scale*Y0[i];// + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*Z0[i]);
//						XYZ[2*NvnG+i] = scale*Z0[i];// + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*Y0[i]);
						XYZ[       i] = scale*(X0[i] + A*dxyz*sin(n*PI/L0*Y0[i])*sin(n*PI/L0*(Z0[i]+0.5)));
						XYZ[1*NvnG+i] = scale*(Y0[i] + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*(Z0[i]+0.5)));
						XYZ[2*NvnG+i] = scale*(Z0[i] + A*dxyz*sin(n*PI/L0*X0[i])*sin(n*PI/L0*Y0[i]));
					}
				} else {
					printf("Error: PeriodicVortex TestCase not supported for dimension d = %d.\n",d), EXIT_MSG;
				}
			} else if (strstr(Geometry,"n-Cylinder")) {
				ToBeCurved_square_to_circle(NvnG,XYZ_S,XYZ);
			} else if (strstr(Geometry,"n-CubeCurved")) {
				if (strstr(DB.GeomSpecifier,"ParabolicPipe"))
					ToBeCurved_parabolic_pipe(NvnG,XYZ_S,XYZ);
				else if (strstr(DB.GeomSpecifier,"EllipticPipe"))
					ToBeCurved_elliptic_pipe(NvnG,XYZ_S,XYZ);
				else if (strstr(DB.GeomSpecifier,"SinusoidalPipe"))
					ToBeCurved_sinusoidal_pipe(NvnG,XYZ_S,XYZ);
				else
					EXIT_UNSUPPORTED;
			}  else {
				printf("%s\n",Geometry);
				printf("Error: Unsupported TestCase for the ToBeCurved MeshType.\n"), EXIT_MSG;
			}

			// Correct internal VOLUME geometry coordinates using blending
			if (nG == 1 && correctTBC && VOLUME->Eclass != C_WEDGE && VOLUME->Eclass != C_PYR)
				correct_ToBeCurved(VOLUME);
		}
	}
}

static void ToBeCurved_square_to_circle(unsigned int Nn, double *XYZ_S, double *XYZ)
{
	/*
	 *	Comments:
	 *		Uses exact r-theta parametrization to transform from square to circular domain.
	 */

	// Initialize DB Parameters
	unsigned int d = DB.d;

	// Standard datatypes
	unsigned int n, dim;
	double       PIo4, XYZn[DMAX], r, t;

	PIo4 = 0.25*PI;
	for (n = 0; n < Nn; n++) {
		for (dim = 0; dim < d; dim++)
			XYZn[dim] = XYZ_S[Nn*dim+n];

		r = array_norm_d(d,XYZn,"Inf");
		t = atan2(XYZn[1],XYZn[0]);

		if      (t >= -1.0*PIo4 && t <  1.0*PIo4) t =          XYZn[1]/r*PIo4;
		else if (t >=  1.0*PIo4 && t <  3.0*PIo4) t = 0.5*PI - XYZn[0]/r*PIo4;
		else if ((t >=  3.0*PIo4 && t <= 4.0*PIo4)|| (t >= -4.0*PIo4 && t < -3.0*PIo4))
		                                          t =     PI - XYZn[1]/r*PIo4;
		else if (t >= -3.0*PIo4 && t < -1.0*PIo4) t = 1.5*PI + XYZn[0]/r*PIo4;
		else
			printf("Error: % .3e % .3e % .3e\n",XYZn[0],XYZn[1],t), EXIT_MSG;

		XYZ[     n] = r*cos(t);
		XYZ[Nn*1+n] = r*sin(t);
		if (d == DMAX)
			XYZ[Nn*2+n] = XYZn[2];
	}
}

static void ToBeCurved_cube_to_sphere(unsigned int Nn, double *XYZ_S, double *XYZ)
{
	/*
	 *	Comments:
	 *		Uses exact parametrization (through polar coordinate transformation) for 2D and Rosca(2011)'s formula for
	 *		the 3D spherical projection.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int dim, n,
	             OrderOut[3];
	int          SignOut;
	double       XYZn[3], XYn[2], XYZn_normInf, beta, *XYZ_Sphere;

	// silence
	for (dim = 0; dim < 3; dim++) {
		if (dim < 2)
			XYn[dim] = 0.0;
		OrderOut[dim] = 0;
	}

	if (d == 3) {
//	if (1||d == 3) {
		for (n = 0; n < Nn; n++) {
			for (dim = 0; dim < d; dim++)
				XYZn[dim] = XYZ_S[Nn*dim+n];

			if (d == 2 || strstr(TestCase,"SupersonicVortex"))
				XYZn[2] = 0.0;

			for (dim = 0; dim < d; dim++) {
				beta = array_norm_d(1,&XYZn[dim],"Inf");
				XYZn_normInf = array_norm_d(3,XYZn,"Inf");
				if (array_norm_diff_d(1,&XYZn_normInf,&beta,"Inf") < EPS)
					break;
			}

			SignOut = sign(XYZn[dim]);
			if (dim == 0) {
				XYn[0] = XYZn[1];
				XYn[1] = XYZn[2];
				OrderOut[0] = 2; OrderOut[1] = 0; OrderOut[2] = 1;
			} else if (dim == 1) {
				XYn[0] = XYZn[0];
				XYn[1] = XYZn[2];
				OrderOut[0] = 0; OrderOut[1] = 2; OrderOut[2] = 1;
			} else if (dim == 2) {
				if (d == 2)
					printf("Error: Invalid entry for d = 2.\n"), EXIT_MSG;

				XYn[0] = XYZn[0];
				XYn[1] = XYZn[1];
				OrderOut[0] = 0; OrderOut[1] = 1; OrderOut[2] = 2;
			}

			XYZ_Sphere = cube_to_sphere(XYn,OrderOut,SignOut,beta); // free

			for (dim = 0; dim < 2; dim++)
				XYZ[Nn*dim+n] = XYZ_Sphere[dim];

			if (d == 3) {
				if (strstr(TestCase,"SupersonicVortex") == NULL)
					XYZ[Nn*2+n] = XYZ_Sphere[2];
				else
					XYZ[Nn*2+n] = XYZ_S[Nn*2+n];
			}

			free(XYZ_Sphere);
		}
	} else if (d == 2) { // Exact r-theta parametrization
		ToBeCurved_square_to_circle(Nn,XYZ_S,XYZ);
	} else {
		printf("Error: Unsupported d.\n"), EXIT_MSG;
	}
}

static double *cube_to_sphere(double XY[2], unsigned int OrderOut[3], int SignOut, double beta)
{
	// Standard datatypes
	unsigned int i;
	double X, Y, XYZ_Sphere_std[3], *XYZ_Sphere;

	X = XY[0];
	Y = XY[1];

	if (array_norm_d(1,&X,"Inf") < 10*EPS && array_norm_d(1,&Y,"Inf") < 10*EPS) {
		for (i = 0; i < 2; i++)
			XY[i] = 0.0;
	} else if (fabs(Y) <= fabs(X)) {
		XY[0] = pow(2.0,0.25)*X/beta*(sqrt(2.0)*cos(Y*PI/(12.0*X))-1.0)/sqrt(sqrt(2.0)-cos(Y*PI/(12.0*X)));
		XY[1] = pow(2.0,0.25)*X/beta*(sqrt(2.0)*sin(Y*PI/(12.0*X)))/sqrt(sqrt(2.0)-cos(Y*PI/(12.0*X)));
	} else if (fabs(X) < fabs(Y)) {
		XY[0] = pow(2.0,0.25)*Y/beta*(sqrt(2.0)*sin(X*PI/(12.0*Y)))/sqrt(sqrt(2.0)-cos(X*PI/(12.0*Y)));
		XY[1] = pow(2.0,0.25)*Y/beta*(sqrt(2.0)*cos(X*PI/(12.0*Y))-1.0)/sqrt(sqrt(2.0)-cos(X*PI/(12.0*Y)));
	} else {
		printf("Error: Invalid condition.\n"), EXIT_MSG;
	}

	X = XY[0];
	Y = XY[1];

	XYZ_Sphere_std[0] = sqrt(1.0-0.25*(pow(X,2.0)+pow(Y,2.0)))*X;
	XYZ_Sphere_std[1] = sqrt(1.0-0.25*(pow(X,2.0)+pow(Y,2.0)))*Y;
	XYZ_Sphere_std[2] = SignOut*(1.0-0.5*(pow(X,2.0)+pow(Y,2.0)));

	XYZ_Sphere = malloc(3 * sizeof *XYZ_Sphere); // keep (requires external free)

	for (i = 0; i < 3; i++)
		XYZ_Sphere[i] = beta*XYZ_Sphere_std[OrderOut[i]];

	return XYZ_Sphere;
}

static void ToBeCurved_TP(unsigned int Nn, double *XYZ_S, double *XYZ)
{
	/*
	 *	Purpose:
	 *		Use Boolean sum projection to linearly blend the effect of curved FACEs into the domain.
	 *
	 *	Comments:
	 *		XZnf and XZnf_match are in row major orientation.
	 *		It is assumed that only one side of element boundaries is curved. This can easily be modified, but would
	 *		require additional arc length computations to locate node positions on other curved sides.
	 *		The convention used below for the vertex and edge numbering was chosen arbitrarily and is not related to
	 *		specifications in 'setup_mesh'.
	 */

	// Initialize DB Parameters
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, dim, n, fn, count, Nv, Ne, Nf,
	             fnCount, NUpdate, ComputeArcLen,
	             xLoc, zLoc, DOrder[2];
	double       XL, XR, ZL, ZR, fNewton, ArcLen, ArcLenTotal,
	             H, s, t, u, sL, sR, tL, tR, uL, uR, sDenom, tDenom, uDenom,
	             XYZn[3], XZnew[2], XYZnew[3], *XZnf, *XZnf_match, *xBounds, *zBounds, *dfdxz_XZnew, *f_XZ, *abcP,
	             F6[2], *Fv, *Fe, *Ff;

	// silence
	XYZn[0] = XYZn[1] = XYZn[2] = 0.0;
	XL = XR = ZL = ZR = 0.0;
	Nv = 0; Ne = 0; Nf = 0;

	XZnf       = malloc (Nn*2 * sizeof *XZnf);       // free
	XZnf_match = malloc (Nn*2 * sizeof *XZnf_match); // free
	xBounds    = malloc(5 * sizeof *xBounds); // free
	zBounds    = malloc(5 * sizeof *zBounds); // free

	if      (d == 2) Nv = 4, Ne = 4,  Nf = 0;
	else if (d == 3) Nv = 8, Ne = 12, Nf = 6;

	Fv = malloc(Nv*d * sizeof *Fv); // free
	Fe = malloc(Ne*d * sizeof *Fe); // free
	Ff = malloc(Nf*d * sizeof *Ff); // free

	fnCount = 0;
	for (n = 0; n < Nn; n++) {
		for (dim = 0; dim < d; dim++)
			XYZn[dim] = XYZ_S[Nn*dim+n];

		if (d == 2)
			XYZn[2] = 0.0;

		if (is_in_blend_region(XYZn)) {
			// Find new coordinates using a proportional arc length parametrization
			XZnew[0] = XYZn[0];
			XZnew[1] = XYZn[2];

			NUpdate = 1;
			if (d == 3)
				NUpdate = 2; // Slower as Nupdate is increased

			ComputeArcLen = 1;
			for (fn = 0; fn < fnCount; fn++) {
				if (array_norm_diff_d(2,&XZnf[fn*2],XZnew,"Inf") < NODETOL) {
					ComputeArcLen = 0;
					for (i = 0; i < 2; i++)
						XZnew[i] = XZnf_match[fn*2+i];
				}
			}

//printf("node, computeArcLen: %d %d\n\n",n,ComputeArcLen);
			if (ComputeArcLen) {
				while (NUpdate--) {
					// Update x-coordinate
					get_blend_bounds(XZnew[0],0.0,&xLoc,&zLoc,xBounds,zBounds);
					XL = xBounds[xLoc];
					XR = xBounds[xLoc+1];
					DOrder[0] = 1; DOrder[1] = 0;
					ArcLenTotal = get_arc_length(XL,XR,XZnew[1],DOrder);

					// Use Newton's method to find updated X
					fNewton = 1.0;
					count = 0;
					while(fabs(fNewton) >= 1e3*EPS) {
						count++;

						ArcLen = get_arc_length(XL,XZnew[0],XZnew[1],DOrder);
						fNewton = ArcLen/ArcLenTotal-(XYZn[0]-XL)/(XR-XL);

						dfdxz_XZnew = eval_TP_function(1,XZnew,DOrder,0,&abcP); // free
						free(abcP);

						XZnew[0] -= ArcLenTotal*fNewton/(sqrt(1+pow(dfdxz_XZnew[0],2.0)));
//printf("%d %d %e %e %e %e\n",n,count,XZnew[0],XL,XR,ArcLen);

						free(dfdxz_XZnew);

						if (count == 1e2)
							printf("Error: Newton's method not converging.\n"), EXIT_MSG;
					}

					if (d == 3) {
						// Update z-coordinate
						get_blend_bounds(0.0,XZnew[1],&xLoc,&zLoc,xBounds,zBounds);
						ZL = zBounds[zLoc];
						ZR = zBounds[zLoc+1];
						DOrder[0] = 0; DOrder[1] = 1;
						ArcLenTotal = get_arc_length(ZL,ZR,XZnew[0],DOrder);

						// Use Newton's method to find updated Z
						fNewton = 1.0;
						count = 0;
						while(fabs(fNewton) >= 1e3*EPS) {
							count++;

							ArcLen = get_arc_length(ZL,XZnew[1],XZnew[0],DOrder);
							fNewton = ArcLen/ArcLenTotal-(XYZn[2]-ZL)/(ZR-ZL);

							dfdxz_XZnew = eval_TP_function(1,XZnew,DOrder,0,&abcP); // free
							free(abcP);

							XZnew[1] -= ArcLenTotal*fNewton/(sqrt(1+pow(dfdxz_XZnew[0],2.0)));
//printf("%d %d %e %e %e %e\n",n,count,XZnew[1],ZL,ZR,ArcLen);

							free(dfdxz_XZnew);

							if (count == 1e2)
								printf("Error: Newton's method not converging.\n"), EXIT_MSG;
						}
					}
				}

				XZnf[fnCount*2+0] = XYZn[0];
				XZnf[fnCount*2+1] = XYZn[2];
				for (i = 0; i < 2; i++)
					XZnf_match[fnCount*2+i] = XZnew[i];

				fnCount++;
			}
//printf("node: %d\n",n);
//array_print_d(1,2,XZnew,'R');

			// Use the Boolean Sum Projection to Find the Final Node Position

			H = 4.0;
			s = XYZn[0]; t = XYZn[1]; u = XYZn[2];

			DOrder[0] = 0; DOrder[1] = 0;
			if (d == 2) {
				f_XZ = eval_TP_function(1,XZnew,DOrder,0,&abcP); // free
				free(abcP);

				sL = XL;  sR = XR;
				tL = 0.0; tR = H;

				Fv[0*d+0] = sL; Fv[0*d+1] = tL;
				Fv[1*d+0] = sR; Fv[1*d+1] = tL;
				Fv[2*d+0] = sR; Fv[2*d+1] = tR;
				Fv[3*d+0] = sL; Fv[3*d+1] = tR;

				Fe[0*d+0] = sL;       Fe[0*d+1] = t;
				Fe[1*d+0] = sR;       Fe[1*d+1] = t;
				Fe[2*d+0] = XZnew[0]; Fe[2*d+1] = f_XZ[0];
				Fe[3*d+0] = s;        Fe[3*d+1] = tR;
				free(f_XZ);

				sDenom = sR-sL;
				tDenom = tR-tL;
				for (dim = 0; dim < d; dim++) {
					XYZnew[dim] = // Vertex Contributions (4)
						- (sR-s)/sDenom*(tR-t)/tDenom*Fv[0*d+dim]
						- (s-sL)/sDenom*(tR-t)/tDenom*Fv[1*d+dim]
						- (s-sL)/sDenom*(t-tL)/tDenom*Fv[2*d+dim]
						- (sR-s)/sDenom*(t-tL)/tDenom*Fv[3*d+dim]
					              // Edge (Face in 2d) Contributions (4)
						+ (sR-s)/sDenom*Fe[0*d+dim]
						+ (s-sL)/sDenom*Fe[1*d+dim]
						+ (tR-t)/tDenom*Fe[2*d+dim]
						+ (t-tL)/tDenom*Fe[3*d+dim];
				}
			} else if (d == 3) {
				f_XZ = eval_TP_function(1,XZnew,DOrder,0,&abcP); // free
				free(abcP);

				sL = XL;  sR = XR;
				tL = 0.0; tR = H;
				uL = ZL;  uR = ZR;

				Fv[0*d+0] = sL; Fv[0*d+1] = tL; Fv[0*d+2] = uL;
				Fv[1*d+0] = sR; Fv[1*d+1] = tL; Fv[1*d+2] = uL;
				Fv[2*d+0] = sL; Fv[2*d+1] = tR; Fv[2*d+2] = uL;
				Fv[3*d+0] = sR; Fv[3*d+1] = tR; Fv[3*d+2] = uL;
				Fv[0*d+0] = sL; Fv[0*d+1] = tL; Fv[4*d+2] = uR;
				Fv[1*d+0] = sR; Fv[1*d+1] = tL; Fv[5*d+2] = uR;
				Fv[2*d+0] = sL; Fv[2*d+1] = tR; Fv[6*d+2] = uR;
				Fv[3*d+0] = sR; Fv[3*d+1] = tR; Fv[7*d+2] = uR;

				Fe[0*d+0]  = s;        Fe[0*d+1]  = tL;      Fe[0*d+2]  = uL;
				Fe[1*d+0]  = s;        Fe[1*d+1]  = tR;      Fe[1*d+2]  = uL;
				Fe[2*d+0]  = XZnew[0]; Fe[2*d+1]  = f_XZ[0]; Fe[2*d+2]  = uR;
				Fe[3*d+0]  = s;        Fe[3*d+1]  = tR;      Fe[3*d+2]  = uR;
				Fe[4*d+0]  = sL;       Fe[4*d+1]  = t;       Fe[4*d+2]  = uL;
				Fe[5*d+0]  = sR;       Fe[5*d+1]  = t;       Fe[5*d+2]  = uL;
				Fe[6*d+0]  = sL;       Fe[6*d+1]  = t;       Fe[6*d+2]  = uR;
				Fe[7*d+0]  = sR;       Fe[7*d+1]  = t;       Fe[7*d+2]  = uR;
				Fe[8*d+0]  = sL;       Fe[8*d+1]  = tL;      Fe[8*d+2]  = u;
				Fe[9*d+0]  = sR;       Fe[9*d+1]  = tL;      Fe[9*d+2]  = u;
				Fe[10*d+0] = sL;       Fe[10*d+1] = tR;      Fe[10*d+2] = u;
				Fe[11*d+0] = sR;       Fe[11*d+1] = tR;      Fe[11*d+2] = u;

				sDenom = sR-sL;
				tDenom = tR-tL;
				for (dim = 0; dim < 2; dim++) {
					F6[dim] = // Vertex Contributions (4)
						- (sR-s)/sDenom*(tR-t)/tDenom*Fv[4*d+dim]
						- (s-sL)/sDenom*(tR-t)/tDenom*Fv[5*d+dim]
						- (sR-s)/sDenom*(t-tL)/tDenom*Fv[6*d+dim]
						- (s-sL)/sDenom*(t-tL)/tDenom*Fv[7*d+dim]
					          // Edge (Face in 2d) Contributions (4)
						+ (tR-t)/tDenom*Fe[2*d+dim]
						+ (t-tL)/tDenom*Fe[3*d+dim]
						+ (sR-s)/sDenom*Fe[6*d+dim]
						+ (s-sL)/sDenom*Fe[7*d+dim];
				}

				Ff[0*d+0] = sL;    Ff[0*d+1] = t;       Ff[0*d+2] = u;
				Ff[1*d+0] = sR;    Ff[1*d+1] = t;       Ff[1*d+2] = u;
				Ff[2*d+0] = s;     Ff[2*d+1] = f_XZ[0]; Ff[2*d+2] = u;
				Ff[3*d+0] = s;     Ff[3*d+1] = tR;      Ff[3*d+2] = u;
				Ff[4*d+0] = s;     Ff[4*d+1] = t;       Ff[4*d+2] = uL;
				Ff[5*d+0] = F6[0]; Ff[5*d+1] = F6[1];   Ff[5*d+2] = uR;

				uDenom = uR-uL;
				for (dim = 0; dim < d; dim++) {
					XYZnew[dim] = // Vertex Contributions (8)
						+ (sR-s)/sDenom*(tR-t)/tDenom*(uR-u)/uDenom*Fv[0*d+dim]
						+ (s-sL)/sDenom*(tR-t)/tDenom*(uR-u)/uDenom*Fv[1*d+dim]
						+ (sR-s)/sDenom*(t-tL)/tDenom*(uR-u)/uDenom*Fv[2*d+dim]
						+ (s-sL)/sDenom*(t-tL)/tDenom*(uR-u)/uDenom*Fv[3*d+dim]
						+ (sR-s)/sDenom*(tR-t)/tDenom*(u-uL)/uDenom*Fv[4*d+dim]
						+ (s-sL)/sDenom*(tR-t)/tDenom*(u-uL)/uDenom*Fv[5*d+dim]
						+ (sR-s)/sDenom*(t-tL)/tDenom*(u-uL)/uDenom*Fv[6*d+dim]
						+ (s-sL)/sDenom*(t-tL)/tDenom*(u-uL)/uDenom*Fv[7*d+dim]
					              // Edge Contributions (12)
						- (tR-t)/tDenom*(uR-u)/uDenom*Fe[0*d+dim]
						- (t-tL)/tDenom*(uR-u)/uDenom*Fe[1*d+dim]
						- (tR-t)/tDenom*(u-uL)/uDenom*Fe[2*d+dim]
						- (t-tL)/tDenom*(u-uL)/uDenom*Fe[3*d+dim]
						- (sR-s)/sDenom*(uR-u)/uDenom*Fe[4*d+dim]
						- (s-sL)/sDenom*(uR-u)/uDenom*Fe[5*d+dim]
						- (sR-s)/sDenom*(u-uL)/uDenom*Fe[6*d+dim]
						- (s-sL)/sDenom*(u-uL)/uDenom*Fe[7*d+dim]
						- (sR-s)/sDenom*(tR-t)/tDenom*Fe[8*d+dim]
						- (s-sL)/sDenom*(tR-t)/tDenom*Fe[9*d+dim]
						- (sR-s)/sDenom*(t-tL)/tDenom*Fe[10*d+dim]
						- (s-sL)/sDenom*(t-tL)/tDenom*Fe[11*d+dim]
					              // Face Contributions (6)
						+ (sR-s)/sDenom*Ff[0*d+dim]
						+ (s-sL)/sDenom*Ff[1*d+dim]
						+ (tR-t)/tDenom*Ff[2*d+dim]
						+ (t-tL)/tDenom*Ff[3*d+dim]
						+ (uR-u)/uDenom*Ff[4*d+dim]
						+ (u-uL)/uDenom*Ff[5*d+dim];
				}
				free(f_XZ);
			}
			for (dim = 0; dim < d; dim++)
				XYZ[Nn*dim+n] = XYZnew[dim];
		} else {
			for (dim = 0; dim < d; dim++)
				XYZ[Nn*dim+n] = XYZn[dim];
		}
	}

	free(XZnf);
	free(XZnf_match);
	free(xBounds);
	free(zBounds);
	free(Fv);
	free(Fe);
	free(Ff);
}

static unsigned int is_in_blend_region(double XYZn[3])
{
	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int d          = DB.d,
	             *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int u1, Inside,
	             IndBound,
	             xLoc, zLoc;
	double       *xBounds, *zBounds, xR, zR;

	u1 = 1;

	xBounds = malloc(5 * sizeof *xBounds); // free
	zBounds = malloc(5 * sizeof *zBounds); // free

	Inside = 0;
	if (strstr(TestCase,"GaussianBump")) {
		get_blend_bounds(XYZn[0],XYZn[2],&xLoc,&zLoc,xBounds,zBounds);
		xR = xBounds[3];
		zR = zBounds[3];

		if (fabs(XYZn[0]) < xR-NODETOL && fabs(XYZn[2]) < zR-NODETOL)
			Inside = u1;
	} else if (strstr(TestCase,"PolynomialBump")) {
		get_blend_bounds(XYZn[0],XYZn[2],&xLoc,&zLoc,xBounds,zBounds);

		if (BumpOrder[0] <= 2)
			IndBound = 4;
		else
			IndBound = 3;

		xR = xBounds[IndBound];
		zR = zBounds[IndBound];

		if (fabs(XYZn[0]) < xR-NODETOL) {
			if (d != 3) {
				Inside = u1;
			} else {
				if (fabs(XYZn[2]) < zR-NODETOL)
					Inside = u1;
			}
		}
	}
	free(xBounds);
	free(zBounds);

	return Inside;
}

static void get_blend_bounds(const double Xn, const double Zn, unsigned int *xLoc, unsigned int *zLoc, double *xBounds,
                             double *zBounds)
{
	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int xLocOut, zLocOut;

	if (strstr(TestCase,"GaussianBump")) {
		xBounds[0] = -100.0; xBounds[1] = -4.0; xBounds[2] = 0.0; xBounds[3] = 4.0; xBounds[4] = 100.0;
		zBounds[0] = -100.0; zBounds[1] = -4.0; zBounds[2] = 0.0; zBounds[3] = 4.0; zBounds[4] = 100.0;
	} else if (strstr(TestCase,"PolynomialBump")) {
		if (BumpOrder[0] <= 2) {
			xBounds[0] = -2.0; xBounds[1] = -1.0; xBounds[2] = 0.0; xBounds[3] = 1.0; xBounds[4] = 2.0;
			zBounds[0] = -2.0; zBounds[1] = -1.0; zBounds[2] = 0.0; zBounds[3] = 1.0; zBounds[4] = 2.0;
		} else {
			xBounds[0] = -100.0; xBounds[1] = -2.0; xBounds[2] = 0.0; xBounds[3] = 2.0; xBounds[4] = 100.0;
			zBounds[0] = -100.0; zBounds[1] = -2.0; zBounds[2] = 0.0; zBounds[3] = 2.0; zBounds[4] = 100.0;
		}
	}

	xLocOut = 0;
	while (Xn > xBounds[xLocOut+1] && xLocOut < 3)
		xLocOut++;

	zLocOut = 0;
	while (Zn > zBounds[zLocOut+1] && zLocOut < 3)
		zLocOut++;

	*xLoc = xLocOut;
	*zLoc = zLocOut;
}

static double get_arc_length(const double XL, const double XR, const double Z, const unsigned int DOrder[2])
{
	/*
	 *	Purpose:
	 *		Compute approximate arc length of supported functions using Romberg's method or exact arc length for
	 *		polynomial functions.
	 *
	 *	Comments:
	 *		Analytical Arc Length Computation:
	 *
	 *			Obtain coefficients for the polynomial under the sqrt for the arc length computation.
	 *			Goal: Compute arc length (AL) of the tensor-product function for fixed x or z (z used in example below)
	 *
	 *			1) f(x,Z) == f(x)*f(Z)
	 *			2) AL =  Integral[sqrt(1+(d/dx{f(x,Z)})^2),{x,XL,XR}]
	 *			      =  Integral[sqrt(1+f(Z)^2*(d/dx{f(x)})^2),{x,XL,XR}]
	 *			      == Integral[sqrt(P(x)),{x,XL,XR}]
	 *			3) Write P(x) =  a*x^2 + b*x + c
	 *			              == a*((x+b1)^2 + c1) if (a != 0)
	 *			4) AL = sqrt(a)/2*((X+b1)*sqrt(P(X)/a) + c1*ln(abs((X+b1)+sqrt(P(X)/a)))) |_{XL,XR}
	 *
	 *			Note: If the order of f(x) is 1, then the integrand is constant and
	 *				AL = sqrt(1+f(Z)^2*(d/dx{f(x)})^2)*(XR-XL)
	 *			Note: If the order of f(x) is greater than 2, then it may be possible to find an anylytical expression
	 *			      for the arc length using elliptic integrals but this has not been investigated.
	 *
	 *	References:
	 *		Wikipedia: https://en.wikipedia.org/wiki/Romberg%27s_method
	 */

	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int PGlobal    = DB.PGlobal,
	             *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int i, AnalyticalArcLen, DOrder_Analytic[2] = {0, 0};
	double       ArcLenOut, XZ[2], *f_XZ, *abcP;

	// silence
	ArcLenOut = 0;
	f_XZ = abcP = NULL;

	AnalyticalArcLen = 0;
	if (strstr(TestCase,"PolynomialBump")) {
		if (array_norm_ui(2,BumpOrder,"Inf") <= 2)
			AnalyticalArcLen = 1;
	}

	if (!AnalyticalArcLen) {
		unsigned int n, l, iMax, lMax, IndX, IndZ,
		             FoundArcLen, Nn;
		double       h, xL, xR,
		             *ArcLen, *XZInt, *dFdX_X, *rst, *w;

		struct S_CUBATURE *CUBDATA = malloc(sizeof *CUBDATA); // free

		set_cubdata(CUBDATA,true,false,"GL",1,PGlobal,cubature_TP); // free

		Nn    = CUBDATA->Nn;
		rst   = CUBDATA->rst;
		w     = CUBDATA->w;

		lMax = 20;

		ArcLen = calloc(pow(lMax,2.0) , sizeof *ArcLen); // free
		XZInt  = calloc(Nn*3          , sizeof *XZInt);  // free

		h = XR-XL;

		if      (DOrder[0] == 1) IndX = 0, IndZ = 1;
		else if (DOrder[1] == 1) IndX = 1, IndZ = 0;
		else                     printf("Error: One entry of DOrder must be 1.\n"), EXIT_MSG;

		for (n = 0; n < Nn; n++) {
			XZInt[Nn*IndZ+n] = Z;
		}

		FoundArcLen = 0;
		for (l = 0; l < lMax && !FoundArcLen; l++) {
			// Compute arc length in each subsection and sum the result
			for (i = 0, iMax = pow(2,l); i < iMax; i++) {
				xL = XL + h*i;
				xR = xL + h;

				for (n = 0; n < Nn; n++)
					XZInt[Nn*IndX+n] = (xL+xR)/2.0 + h/2.0*rst[n];

				dFdX_X = eval_TP_function(Nn,XZInt,DOrder,0,&abcP); // free
				free(abcP);

				for (n = 0; n < Nn; n++)
					ArcLen[l*lMax] += h/2.0*w[n]*sqrt(1.0+pow(dFdX_X[n],2.0));

				free(dFdX_X);
			}

			// Compute Richardson Extrapolation
			for (i = 1; i <= l; i++) {
				ArcLen[l*lMax+i] =
					(pow(4.0,(double) i)*ArcLen[l*lMax+i-1]-ArcLen[(l-1)*lMax+i-1])/(pow(4.0,(double) i)-1.0);

				if (l > 3 && array_norm_diff_d(1,&ArcLen[l*lMax+i],&ArcLen[l*lMax+i-1],"Inf") < NODETOL) {
					ArcLenOut = ArcLen[l*lMax+i];
					FoundArcLen = 1;
					break;
				}
			}
			h /= 2.0;
		}
		if (!FoundArcLen)
			printf("Error: Arc length was not found using lMax = %d.\n",lMax), EXIT_MSG;
//array_print_d(lMax,lMax,ArcLen,'R');

		free(ArcLen);
		free(XZInt);
		free(rst);
		free(w);
		free(CUBDATA);
	} else {
		double aF, bF, a, b, c, b1, c1, X, P_X;

		// Obtain polynomial coefficients of varying portion of f(x,z)
		if (DOrder[0] == 1) {
			XZ[0] = XR; XZ[1] = Z;
			f_XZ = eval_TP_function(1,XZ,DOrder_Analytic,1,&abcP); // free
		} else if (DOrder[1] == 1) {
			XZ[0] = Z; XZ[1] = XR;
			f_XZ = eval_TP_function(1,XZ,DOrder_Analytic,3,&abcP); // free
		}
		free(f_XZ);

		// Find coefficients of P(x)
		// Note: As the functions are symmetric about (0,0), some simplifications have been implemented below which take
		//       this into account.

		aF = abcP[0]; bF = abcP[1];
		free(abcP);

		a = pow(2.0*aF,2.0); b = 4.0*aF*bF; c = 1.0+pow(bF,2.0);
		b1 = b/(2.0*a); c1 = c/a-pow(b1,2.0);

		// Compute Arc Length using the analytical expression
		if (fabs(aF) < EPS) {
			ArcLenOut = sqrt(c)*(XR-XL);
		} else {
			ArcLenOut = 0.0;
			for (i = 0; i < 2; i++) {
				if (i == 0) X = XL;
				else        X = XR;

				P_X = a*pow(X,2.0) + b*X + c;

				ArcLenOut += pow(-1.0,i+1.0)*sqrt(a)/2.0*((X+b1)*sqrt(P_X/a)+c1*log(fabs((X+b1)+sqrt(P_X/a))));
			}
		}
	}

	return ArcLenOut;
}

static double *eval_TP_function(const unsigned int Nn, const double *XZ, const unsigned int DOrder[2],
                                const unsigned int Single, double **abcP)
{
	/*
	 *	Purpose:
	 *		Evaluate supported (T)ensor-(P)roduct functions.
	 */

	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int n, dim,
	             xLoc, zLoc;
	double       X, Z, F_X, F_Z, a, b, c, e, h,
	             abcX[3], abcZ[3], *Output, *xBounds, *zBounds, *abcOut;

	// silence
	a = b = c = 0;
	abcX[0] = abcX[1] = abcX[2] = 0;
	abcZ[0] = abcZ[1] = abcZ[2] = 0;
	F_X = F_Z = 0;

	Output  = malloc(Nn * sizeof *Output);  // keep (requires external free)
	xBounds = malloc(5  * sizeof *xBounds); // free
	zBounds = malloc(5  * sizeof *zBounds); // free
	abcOut  = calloc(3  , sizeof *abcOut);  // keep (requires external free)

	for (n = 0; n < Nn; n++) {
		X = XZ[0*Nn+n];
		Z = XZ[1*Nn+n];

		get_blend_bounds(X,Z,&xLoc,&zLoc,xBounds,zBounds);

		// Compute function output in the correct region
		if (strstr(TestCase,"GaussianBump")) {
			// f2(X) = f3(X) = a*exp(-(X-b).^2/(2*c^2))
			// f2(Z) = f3(Z) = a*exp(-(Z-b).^2/(2*c^2))
			// f1(X),f4(X),f1(Z),f4(Z): Not Needed

			a = 0.5;
			b = 0.0;
			c = 0.6;

			if (xLoc == 1 || xLoc == 2) {
				if      (DOrder[0] == 0) F_X = a*exp(-1.0*pow(X-b,2.0)/(2.0*pow(c,2.0)));
				else if (DOrder[0] == 1) F_X = a*(b-X)/pow(c,2.0)*exp(-1.0*pow(X-b,2.0)/(2.0*pow(c,2.0)));
				else
					printf("Error: Add in support for higher derivatives.\n"), EXIT_MSG;
			} else {
				printf("Error: Invalid function evaluation for GaussianBump.\n"), EXIT_MSG;
			}

			if (zLoc == 1 || zLoc == 2) {
				if      (DOrder[1] == 0) F_Z = a*exp(-1.0*pow(Z-b,2.0)/(2.0*pow(c,2.0)));
				else if (DOrder[1] == 1) F_Z = a*(b-Z)/pow(c,2.0)*exp(-1.0*pow(Z-b,2.0)/(2.0*pow(c,2.0)));
				else
					printf("Error: Add in support for higher derivatives.\n"), EXIT_MSG;
			} else {
				printf("Error: Invalid function evaluation for GaussianBump.\n"), EXIT_MSG;
			}
		} else if (strstr(TestCase,"PolynomialBump")) {
			h = 0.5;

			if (array_norm_ui(2,BumpOrder,"Inf") <= 2) {
				// Bumps in x-direction
				if (xLoc == 0) {
					// Bump 1
					if      (BumpOrder[1] == 0) abcX[0] = 0.0  , abcX[1] = 0.0  , abcX[2] = 0.0  ; // f1(X) = 0
					else if (BumpOrder[1] == 1) abcX[0] = 0.0  , abcX[1] = h/2.0, abcX[2] = h    ; // f1(X) = h/2*(2+X)
					else if (BumpOrder[1] == 2) abcX[0] = h/2.0, abcX[1] = 2.0*h, abcX[2] = 2.0*h; // f1(X) = h/2*(2+X)^2
				} else if (xLoc == 1 || xLoc == 2) {
					// Bumps 2 and 3 (f2(X) == f3(X))
					if (BumpOrder[1] == 0) abcX[0] = -h    , abcX[1] = 0.0, abcX[2] = h; // f2(X) = h*(1-X^2)
					else                   abcX[0] = -0.5*h, abcX[1] = 0.0, abcX[2] = h; // f2(X) = h*(1-X^2/2)
				} else if (xLoc == 3) {
					// Bump 4
					if      (BumpOrder[1] == 0) abcX[0] = 0.0  , abcX[1] = 0.0   , abcX[2] = 0.0  ; // f4(X) = 0
					else if (BumpOrder[1] == 1) abcX[0] = 0.0  , abcX[1] = -h/2.0, abcX[2] = h    ; // f4(X) = h/2*(2-X)
					else if (BumpOrder[1] == 2) abcX[0] = h/2.0, abcX[1] = -2.0*h, abcX[2] = 2.0*h; // f4(X) = h/2*(2-X)^2
				}

				if      (DOrder[0] == 0) a = abcX[0], b = abcX[1]    , c = abcX[2];
				else if (DOrder[0] == 1) a = 0.0    , b = 2.0*abcX[0], c = abcX[1];

				F_X = a*pow(X,2.0) + b*X + c;

				// Bumps in z-direction
				if (zLoc == 0) {
					// Bump 1
					if      (BumpOrder[1] == 0) abcZ[0] = 0.0  , abcZ[1] = 0.0  , abcZ[2] = 0.0  ; // f1(Z) = 0
					else if (BumpOrder[1] == 1) abcZ[0] = 0.0  , abcZ[1] = h/2.0, abcZ[2] = h    ; // f1(Z) = h/2*(2+Z)
					else if (BumpOrder[1] == 2) abcZ[0] = h/2.0, abcZ[1] = 2.0*h, abcZ[2] = 2.0*h; // f1(Z) = h/2*(2+Z)^2
				} else if (zLoc == 1 || zLoc == 2) {
					// Bumps 2 and 3 (f2(Z) == f3(Z))
					if (BumpOrder[1] == 0) abcZ[0] = -h    , abcZ[1] = 0.0, abcZ[2] = h; // f2(Z) = h*(1-Z^2)
					else                   abcZ[0] = -0.5*h, abcZ[1] = 0.0, abcZ[2] = h; // f2(Z) = h*(1-Z^2/2)
				} else if (zLoc == 3) {
					// Bump 4
					if      (BumpOrder[1] == 0) abcZ[0] = 0.0  , abcZ[1] = 0.0   , abcZ[2] = 0.0  ; // f4(Z) = 0
					else if (BumpOrder[1] == 1) abcZ[0] = 0.0  , abcZ[1] = -h/2.0, abcZ[2] = h    ; // f4(Z) = h/2*(2-Z)
					else if (BumpOrder[1] == 2) abcZ[0] = h/2.0, abcZ[1] = -2.0*h, abcZ[2] = 2.0*h; // f4(Z) = h/2*(2-Z)^2
				}

				if      (DOrder[1] == 0) a = abcZ[0], b = abcZ[1]    , c = abcZ[2];
				else if (DOrder[1] == 1) a = 0.0    , b = 2.0*abcZ[0], c = abcZ[1];
				F_Z = a*pow(Z,2.0) + b*Z + c;
			} else {
				// f2(X) = f3(X) = h*(1-1/2*X^2+1/16*X^4)
				// f2(Z) = f3(Z) = h*(1-1/2*Z^2+1/16*Z^4)
				// f1(X), f4(X), f1(Z), f4(Z): Not Needed

				// Note: No odd components
				a = 1;
				c = -0.5;
				e = 0.0625;

				if (xLoc == 1 || xLoc == 2) {
					if      (DOrder[0] == 0) F_X = h*(a+c*pow(X,2.0)+e*pow(X,4.0));
					else if (DOrder[0] == 1) F_X = h*(2.0*c*X+4.0*e*pow(X,3.0));
				} else {
					printf("Error: eval_TP_function should not be needed in this region for the high-order Bump.\n");
					EXIT_MSG;
				}

				if (zLoc == 1 || zLoc == 2) {
					if      (DOrder[1] == 0) F_Z = h*(a+c*pow(Z,2.0)+e*pow(Z,4.0));
					else if (DOrder[1] == 1) F_Z = h*(2.0*c*Z+4.0*e*pow(Z,3.0));
				} else {
					printf("Error: eval_TP_function should not be needed in this region for the high-order Bump.\n");
					EXIT_MSG;
				}
			}
		} else {
			printf("Error: Unsupported TP_function.\n"), EXIT_MSG;
		}
		if (!Single) {
			Output[n] = F_X*F_Z;
		} else {
			if (Single == 1) {
				Output[n] = F_X;
				for (dim = 0; dim < 3; dim++)
					abcOut[dim] = abcX[dim];
			} else if (Single == 3) {
				Output[n] = F_Z;
				for (dim = 0; dim < 3; dim++)
					abcOut[dim] = abcZ[dim];
			} else {
				printf("Error: Invalid value for Single.\n"), EXIT_MSG;
			}
		}
	}
	free(xBounds);
	free(zBounds);

	*abcP = abcOut;
	return Output;
}
