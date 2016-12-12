// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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

typedef void (*compute_pc_tdef) (struct S_pc *data);
typedef void (*compute_XYZ_tdef) (struct S_XYZ *data);

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

	if (norm_rIn < 2e-1) {
		r = rIn;
	} else if (norm_rOut < 2e-1) {
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
	double       r, D, ABC[3], XYZ[DMAX] = {0.0}, XYZ_C[DMAX] ={0.0};

	Indn = 0;
	if (strstr(Geometry,"dm1-Spherical_Section")) {
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
			for (dim = 0; dim < DMAX; dim++) {
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
	} else if (strstr(Geometry,"Ringleb")) {
		if (d != 2)
			printf("Error: Unsupported.\n");

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

				k = 0.5*(KMin+KMax);
				for (i = 0, iMax = 100; i < iMax; i++) {
					// Find point on Ringleb surface for guessed k
					a   = sqrt(1.0-0.5*GM1*q*q);
					rho = pow(a,2.0/GM1);
					J   = 1.0/a+1.0/(3.0*pow(a,3.0))+1.0/(5.0*pow(a,5.0))-0.5*log((1.0+a)/(1.0-a));

					x = 1.0/(2.0*rho)*(2.0/(k*k)-1.0/(q*q))-0.5*J;
					y = sign_y/(k*rho*q)*sqrt(1.0-pow(q/k,2.0));

					// Find the point on the normal line which is closest to the located point on the boundary
					alpha =  nx*(x-xS)+ny*(y-yS);
					beta  = -ny*(x-xS)+nx*(y-yS);

					// Find new point on the surface based on x-coordinate found above and update k
					x = xS + nx*alpha;
					k = sqrt(2.0/(2.0*rho*(x+0.5*J)+1.0/(q*q)));

					if (fabs(beta) < 1e1*EPS)
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
				Ringleb_boundary(&x,&y,0.0,k,'w');
				if (fabs(y-yS) < NODETOL_MESH)
					aStep = 0.0;

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
			printf("Error: Unsupported.\n");
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

			p = acos(Z/rXYZ);
		}

		PComps[0*Nn+n] = t;
		PComps[1*Nn+n] = p;
	}

	if (d == 3) {
		// If one (and only one) of the vertices lies on a pole of the sphere, re-calculate theta for the vertex based
		// on the theta values of the other vertices.
		for (n = 0; n < Nn; n++) {
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
			printf("Error: Unsupported.\n");
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

static void select_functions_Curved(compute_pc_tdef *compute_pc, compute_XYZ_tdef *compute_XYZ)
{
	// Initialize DB Parameters
	char *Geometry = DB.Geometry;

	if (strstr(Geometry,"dm1-Spherical_Section")) {
		*compute_pc  = compute_pc_dsphere;
		*compute_XYZ = compute_XYZ_dsphere;
	} else if (strstr(Geometry,"Ringleb")) {
		// Write these functions if ARC LENGTH projection is desired.
		*compute_pc  = NULL;
		*compute_XYZ = NULL;
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

static double *compute_BlendV(struct S_Blend *data)
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

	if (Blending == HESTHAVEN && type == TRI) {
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
			BlendV[n] = pow(1.0-I_vGs_vGc[n*Nve+b],2.0);
		}
	} else if (Blending == SZABO_BABUSKA || EclassV == C_SI) {
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
	} else if (Blending == GORDON_HALL && EclassV == C_TP) {
		for (n = 0; n < NvnG; n++) {
			BlendV[n] = 0.0;
			for (ve = 0; ve < Nbve[b]; ve++)
				BlendV[n] += I_vGs_vGc[n*Nve+VeBcon[b*NbveMax+ve]];
		}
	} else {
		printf("Error: Unsupported.\n");
	}
/*
if (b == 2) {
printf("%d\n",b);
array_print_d(NvnG,1,BlendV,'R');
array_print_d(NvnG,Nve,I_vGs_vGc,'R');
array_print_d(NvnG,Nbve[b],I_bGs_vGc,'R');
array_print_d(NvnG,data->NbnG,data->I_bGc_vGc,'R');
//EXIT_MSG;
}
*/
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
	I_vGc_bGc = data->I_vGc_bGc;
	I_bGc_vGc = data->I_bGc_vGc;

	data_pc  = malloc(sizeof *data_pc);  // free
	data_XYZ = malloc(sizeof *data_XYZ); // free

	compute_pc  = data->compute_pc;
	compute_XYZ = data->compute_XYZ;

	if (Parametrization == ARC_LENGTH) {
		I_vGs_vGc = data->I_vGs_vGc;

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
	} else if (Parametrization == NORMAL) {
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
		array_free2_d(DMAX,VeXYZ);

		// Compute normal distance from straight FACE to curved geometry
		compute_normal_displacement(NbnG,0,XYZ_S,n_S,XYZ_CmS,BC);

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
	unsigned int d         = DB.d,
	             *PGc      = DB.PGc,
	             Blending  = DB.Blending;

	// Standard datatypes
	unsigned int b, ve, dim, n, Vb, PV, IndSurf, Vtype,
	             Nve, Nb, *Nbve, *VeBcon, NbveMax, NbrefMax, NvnG, BC, *VeInfo;
	double       *XYZ, *XYZ_update, *XYZ_vV, *BlendV, *BlendVe, *I_vGs_vGc;

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
		if (PGc[2] == 2)
			PV = 2;
		else if (PGc[1] == 2)
			PV = 1;
		else
			printf("Error: Add specific blending operators for vertex blending for PGc != P or P+1.\n"), EXIT_MSG;

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
			data_blend->I_vGc_bGc = ELEMENT->I_vGc_eGc[PV][PV][Vb];
			data_blend->I_bGc_vGc = ELEMENT->I_eGc_vGc[PV][PV][Vb];
			data_blend->I_bGs_vGc = ELEMENT->I_eGs_vGc[1][PV][Vb];
		} else if (BType == 'f') {
			ELEMENT_B = get_ELEMENT_F_type(ELEMENT->type,b);
			data_blend->I_vGc_bGc = ELEMENT->I_vGc_fGc[PV][PV][Vb];
			data_blend->I_bGc_vGc = ELEMENT->I_fGc_vGc[PV][PV][Vb];
			data_blend->I_bGs_vGc = ELEMENT->I_fGs_vGc[1][PV][Vb];
		}

		data_blend->b    = b;
		data_blend->NbnG = ELEMENT_B->NvnGc[PV];
		data_blend->BC   = BC;

		BlendV     = compute_BlendV(data_blend);     // free
		XYZ_update = compute_XYZ_update(data_blend); // free

		// Blend BOUNDARY perturbation to VOLUME geometry nodes
		for (n = 0; n < NvnG; n++) {
			if (BlendV[n] < EPS)
				continue;

			for (dim = 0; dim < d; dim++)
				XYZ[n+NvnG*dim] += BlendV[n]*XYZ_update[n+NvnG*dim];
		}
		free(XYZ_update);
		free(BlendV);
	}

	// Add additional contribution from vertices if applicable
	if (0&&Blending == NIELSON && Vtype == TRI) {
		BlendVe = malloc(NvnG*Nve * sizeof *BlendVe); // free

		for (n = 0; n < NvnG; n++) {
		for (ve = 0; ve < Nve; ve++) {
			BlendVe[n*Nve+ve] = pow(I_vGs_vGc[n*Nve+ve],2.0);
		}}
		mm_d(CBCM,CBT,CBNT,NvnG,d,Nve,0.0,1.0,BlendVe,XYZ_vV,XYZ);
		free(BlendVe);
	} else {
		// Do nothing
	}

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
