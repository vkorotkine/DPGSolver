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
 *		Gordon(1973)-Construction_of_Curvilinear_Coordinate_Systems_and_Applications_to_Mesh_Generation
 *		Gordon(1973)-Transfinite_Element_Methods-Blending-Function_Interpolation_over_Arbitraty_Curved_Element_Domains
 *		Rosca(2011)-Uniform_Spherical_Grids_via_Equal_Area_Projection_from_the_Cube_to_the_Sphere
 */

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
		printf("Error: Invalid condition in cube_to_sphere.\n"), exit(1);
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

static void ToBeCurved_cube_to_sphere(unsigned int Nn, double *XYZs, double *XYZ)
{
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

	for (n = 0; n < Nn; n++) {
		for (dim = 0; dim < d; dim++)
			XYZn[dim] = XYZs[Nn*dim+n];

		if (d == 2 || strstr(TestCase,"SupersonicVortex") != NULL)
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
				printf("Error: Invalid entry in ToBeCurved_cube_to_sphere for d = 2.\n"), exit(1);

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
				XYZ[Nn*2+n] = XYZs[Nn*2+n];
		}
		
		free(XYZ_Sphere);
	}
}

static void get_blend_bounds(const double Xn, const double Zn, unsigned int *xLoc, unsigned int *zLoc, double *xBounds,
                             double *zBounds)
{
	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int xLocOut, zLocOut;

	if (strstr(TestCase,"GaussianBump") != NULL) {
		xBounds[0] = -100.0; xBounds[1] = -4.0; xBounds[2] = 0.0; xBounds[3] = 4.0; xBounds[4] = 100.0;
		zBounds[0] = -100.0; zBounds[1] = -4.0; zBounds[2] = 0.0; zBounds[3] = 4.0; zBounds[4] = 100.0;
	} else if (strstr(TestCase,"PolynomialBump") != NULL) {
		if (BumpOrder[0] == 0) {
			xBounds[0] = -100.0; xBounds[1] = -1.0; xBounds[2] = 0.0; xBounds[3] = 1.0; xBounds[4] = 100.0;
			zBounds[0] = -100.0; zBounds[1] = -1.0; zBounds[2] = 0.0; zBounds[3] = 1.0; zBounds[4] = 100.0;
		} else if (BumpOrder[0] <= 2) {
			xBounds[0] = -2.0; xBounds[1] = -1.0; xBounds[2] = 0.0; xBounds[3] = 1.0; xBounds[4] = 2.0;
			zBounds[0] = -2.0; zBounds[1] = -1.0; zBounds[2] = 0.0; zBounds[3] = 1.0; zBounds[4] = 2.0;
		} else {
			xBounds[0] = -100.0; xBounds[1] = -2.0; xBounds[2] = 0.0; xBounds[3] = 2.0; xBounds[4] = 100.0;
			zBounds[0] = -100.0; zBounds[1] = -2.0; zBounds[2] = 0.0; zBounds[3] = 2.0; zBounds[4] = 100.0;
		}
	}

	xLocOut = 0;
	while (Xn > xBounds[xLocOut] && xLocOut < 3)
		xLocOut++;

	zLocOut = 0;
	while (Zn > zBounds[zLocOut] && zLocOut < 3)
		zLocOut++;
	
	*xLoc = xLocOut;
	*zLoc = zLocOut;
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
	if (strstr(TestCase,"GaussianBump") != NULL) {
		Inside = u1;
	} else if (strstr(TestCase,"PolynomialBump") != NULL) {
		get_blend_bounds(XYZn[0],XYZn[2],&xLoc,&zLoc,xBounds,zBounds);

		if (BumpOrder[0] > 0 && BumpOrder[0] < 3)
			IndBound = 4;
		else
			IndBound = 3;

		xR = xBounds[IndBound];
		zR = zBounds[IndBound];

		if (array_norm_diff_d(1,&XYZn[0],&xR,"Inf") < NODETOL) {
			if (d != 3) {
				Inside = u1;
			} else {
				if (array_norm_diff_d(1,&XYZn[2],&zR,"Inf") < NODETOL)
					Inside = u1;
			}
		}
	}
	free(xBounds);
	free(zBounds);

	return Inside;
}

static double *eval_TP_function(const unsigned int Nn, const double *XYZ, const unsigned int DOrder[2],
                                const unsigned int Single)
{
	// Initialize DB Parameters
	char         *TestCase  = DB.TestCase;
	unsigned int *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int n,
	             xLoc, zLoc;
	double       X, Z, F_X, F_Z,
	             *Output, *xBounds, *zBounds;

	Output  = malloc(Nn * sizeof *Output);  // keep (requires external free)
	xBounds = malloc(5  * sizeof *xBounds); // free
	zBounds = malloc(5  * sizeof *zBounds); // free

	for (n = 0; n < Nn; n++) {
		X = XYZ[0*Nn+n];
		Z = XYZ[2*Nn+n];

		get_blend_bounds(X,Z,&xLoc,&zLoc,xBounds,zBounds);

		// Compute function output in the correct region
		if (strstr(TestCase,"GaussianBump") != NULL) {
			// f2(X) = f3(X) = a*exp(-(X-b).^2/(2*c^2))
			// f2(Z) = f3(Z) = a*exp(-(Z-b).^2/(2*c^2))
			// f1(X),f4(X),f1(Z),f4(Z): Not Needed

			double a, b, c;

			a = 0.5;
			b = 0.0;
			c = 0.6;

			if (xLoc == 1 || xLoc == 2) {
				if      (DOrder[0] == 0) F_X = a*exp(-1.0*pow(X-b,2.0)/(2.0*pow(c,2.0)));
				else if (DOrder[0] == 1) F_X = a*(b-X)/pow(c,2.0)*exp(-1.0*pow(X-b,2.0)/(2.0*pow(c,2.0)));
				else
					printf("Error: Add in support for higher derivatives.\n"), exit(1);
			} else {
				printf("Error: Invalid function evaluation for GaussianBump.\n"), exit(1);
			}

			if (zLoc == 1 || zLoc == 2) {
				if      (DOrder[1] == 0) F_Z = a*exp(-1.0*pow(Z-b,2.0)/(2.0*pow(c,2.0)));
				else if (DOrder[1] == 1) F_Z = a*(b-Z)/pow(c,2.0)*exp(-1.0*pow(Z-b,2.0)/(2.0*pow(c,2.0)));
				else
					printf("Error: Add in support for higher derivatives.\n"), exit(1);
			} else {
				printf("Error: Invalid function evaluation for GaussianBump.\n"), exit(1);
			}
		} else if (strstr(TestCase,"PolynomialBump") != NULL) {
			double h;

			h = 0.5;

			printf("Finish eval_TP_function for polynomialbump.\n");
			exit(1);
		} else {
			printf("Error: Unsupported TP_function.\n"), exit(1);
		}
		if (!Single) {
			Output[n] = F_X*F_Z;
		} else {
			if      (Single == 1) Output[n] = F_X;
			else if (Single == 3) Output[n] = F_Z;
			else
				printf("Error: Invalid value for Single in eval_TP_function.\n"), exit(1);
		}
	}
	free(xBounds);
	free(zBounds);

	return Output;
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
	unsigned int P          = DB.P,
	             *BumpOrder = DB.BumpOrder;

	// Standard datatypes
	unsigned int AnalyticalArcLen;
	double       ArcLenOut;

	AnalyticalArcLen = 0;
	if (strstr(TestCase,"PolynomialBump") != NULL) {
		if (array_norm_ui(2,BumpOrder,"Inf") <= 2)
			AnalyticalArcLen = 1;
	}

	if (!AnalyticalArcLen) {
		unsigned int i, n, l, iMax, lMax, IndX, IndZ,
		             FoundArcLen, Nn, Ns, *symms;
		double       h, xL, xR,
		             *ArcLen, *XYZInt, *XInt, *dFdX_X, *rst, *w;

		cubature_TP(&rst,&w,&symms,&Nn,&Ns,1,P,1,"GL"); // free

		lMax = 10;
// lMax = 20

		ArcLen = malloc(pow(lMax,2.0) * sizeof *ArcLen); // free
		XYZInt = calloc(Nn*3          , sizeof *XYZInt); // free
		XInt   = malloc(Nn*1          * sizeof *XInt);   // free

		h = XR-XL;

		if      (DOrder[0] == 1) IndX = 0, IndZ = 2;
		else if (DOrder[1] == 1) IndX = 2, IndZ = 0;
		else                     printf("Error: One entry of DOrder must be 1.\n"), exit(1);

		for (n = 0; n < Nn; n++) {
			XYZInt[Nn*1+n]    = 0.0;
			XYZInt[Nn*IndZ+n] = Z;
		}

		FoundArcLen = 0;
		for (l = 0; l < lMax && !FoundArcLen; l++) {
			// Compute arc length in each subsection and sum the result
			for (i = 0, iMax = pow(2,l); i < iMax; i++) {
				xL = XL + h*i;
				xR = xL + h;

				for (n = 0; n < Nn; n++)
					XYZInt[Nn*IndX+n] = (xL+xR)/2.0 + h/2.0*rst[n];

				dFdX_X = eval_TP_function(Nn,XYZInt,DOrder,0); // free

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

		free(ArcLen);
		free(rst);
		free(w);
		free(symms);
	} else {

	}

	ArcLenOut = 0.0;
	return ArcLenOut;
}

static void ToBeCurved_TP(unsigned int Nn, double *XYZs, double *XYZ)
{
	/*
	 *	Purpose:
	 *		Use Boolean sum projection to linearly blend the effect of curved FACETs into the domain.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, dim, n, fn, update, count,
	             fnCount, NUpdate, ComputeArcLen,
	             xLoc, zLoc, DOrder[2];
	double       XYZn[3], XZnew[2], *XZnf, *XZnf_match, *xBounds, *zBounds,
	             XC, ZC, XL, XR, ZL, ZR, xL, xR, fNewton, ArcLenTotal;

	XZnf       = malloc (Nn*d * sizeof *XZnf);       // free
	XZnf_match = malloc (Nn*d * sizeof *XZnf_match); // free
	xBounds    = malloc(5 * sizeof *xBounds); // free
	zBounds    = malloc(5 * sizeof *zBounds); // free

array_print_d(Nn,d,XYZs,'C');

	fnCount = 0;
	for (n = 0; n < Nn; n++) {
		for (dim = 0; dim < d; dim++)
			XYZn[dim] = XYZs[Nn*dim+n];

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
						XZnew[i] = XZnf[fn*2+i];
				}
			}

			if (ComputeArcLen) {
				while (NUpdate--) {
					// Update x-coordinate
					ZC = XZnew[1];

					get_blend_bounds(XZnew[0],0.0,&xLoc,&zLoc,xBounds,zBounds);
					XL = xBounds[xLoc];
					XR = xBounds[xLoc+1];
					DOrder[0] = 1; DOrder[1] = 0;
					ArcLenTotal = get_arc_length(XL,XR,ZC,DOrder);

					// Use Newton's method to find updated X
					fNewton = 1.0;
					count = 0;

					while(fabs(fNewton) >= 1e3*EPS) {
						count++;


						if (count == 1e2)
							printf("Error: Newton's method not converging in ToBeCurved_TP.\n"), exit(1);
					}

					exit(1);
				}


			}



		}
	}

	free(XZnf);
	free(XZnf_match);
	free(xBounds);
	free(zBounds);
}

void setup_ToBeCurved(struct S_VOLUME *VOLUME)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, dim,
	             NvnG;
	double *XYZ, *XYZs;

	NvnG = VOLUME->NvnG;
	XYZs = VOLUME->XYZs;

	XYZ = malloc (NvnG*d * sizeof *XYZ); // keep
	VOLUME->XYZ = XYZ;

	if (strstr(TestCase,"dSphericalBump")   != NULL ||
	    strstr(TestCase,"PorousdSphere")    != NULL ||
	    strstr(TestCase,"SupersonicVortex") != NULL) {
			ToBeCurved_cube_to_sphere(NvnG,XYZs,XYZ);
	} else if (strstr(TestCase,"GaussianBump")     != NULL ||
	           strstr(TestCase,"PolynomialBump")   != NULL) {
array_print_d(NvnG,d,XYZs,'C');
			ToBeCurved_TP(NvnG,XYZs,XYZ);

array_print_d(NvnG,d,XYZ,'C');
exit(1);

	} else if (strstr(TestCase,"PeriodicVortex") != NULL) {
		double n = 2.0, A = 0.1, L0 = 2.0, dxyz = 1.0;
		double *X0, *Y0, *Z0;

		// silence
		Z0 = NULL;

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
	} else {
		printf("Error: Unsupported TestCase for the ToBeCurved MeshType.\n"), exit(1);
	}
}
