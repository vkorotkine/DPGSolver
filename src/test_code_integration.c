// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_code_integration.h"

#include <string.h>

#include "petscsys.h"
#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "memory_free.h"
#include "array_norm.h"
#include "array_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void update_MeshFile(void)
{
	// Standard datatypes
	char *d, *ML;

	d  = malloc(STRLEN_MIN * sizeof *d);  // free
	ML = malloc(STRLEN_MIN * sizeof *ML); // free

	sprintf(d,"%d",DB.d);
	sprintf(ML,"%d",DB.ML);

	strcpy(DB.MeshFile,"");
	strcat(DB.MeshFile,DB.MeshPath);
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.TestCase);
	strcat(DB.MeshFile,strcat(d,"D_"));
	strcat(DB.MeshFile,DB.MeshType);
	strcat(DB.MeshFile,strcat(ML,"x.msh"));

	free(d);
	free(ML);
}

void code_startup(int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize(&nargc,&argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	initialization(nargc,argv);
	if (update_argv) {
		strcpy(DB.TestCase,TestDB.TestCase);
		DB.PGlobal = TestDB.PGlobal;
		DB.ML      = TestDB.ML;
		update_MeshFile();
	}

	setup_parameters();
	if (update_argv)
		setup_parameters_L2proj();

	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();

	initialize_test_case(Nref);
}

void code_cleanup(void)
{
	mesh_to_level(0);
	memory_free();
}

void check_convergence_orders(const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                              const unsigned int PMax, unsigned int *pass)
{
	// Initialize DB Parameters
	char         *TestCase = TestDB.TestCase,
	             *MeshType = DB.MeshType;
	unsigned int d         = DB.d;

	// Standard datatypes
	char         f_name[STRLEN_MAX], string[STRLEN_MAX], StringRead[STRLEN_MAX], *data;
	unsigned int i, ML, P, NVars, NML, NP, Indh, *VarsToCheck, *OrderIncrement;
	int          offset;
	double       **L2Errors, **ConvOrders, *h, tmp_d;

	FILE *fID;

	// silence
	NVars = 0;

	if (strstr(TestCase,"Poisson")) {
		NVars = DMAX+1;
	} else if (strstr(TestCase,"SupersonicVortex")) {
		NVars = DMAX+2+1;
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}

	VarsToCheck    = malloc(NVars * sizeof *VarsToCheck);    // free
	OrderIncrement = malloc(NVars * sizeof *OrderIncrement); // free
	if (strstr(TestCase,"Poisson")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 0;
			if (i == 0) {
				OrderIncrement[i] = 1;
				VarsToCheck[i]    = 1;
			} else if (i <= d) {
				VarsToCheck[i]    = 1;
			} else {
				VarsToCheck[i]    = 0;
			}
		}
	} else if (strstr(TestCase,"SupersonicVortex")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 1;
			if (i <= d || i > DMAX) {
				VarsToCheck[i] = 1;
			} else {
				VarsToCheck[i] = 0;
			}
		}
	}


	NML = MLMax-MLMin+1;
	NP  = PMax-PMin+1;

	L2Errors   = malloc(NVars * sizeof *L2Errors);   // free
	ConvOrders = malloc(NVars * sizeof *ConvOrders); // free
	for (i = 0; i < NVars; i++) {
		L2Errors[i]   = calloc(NML*NP , sizeof *L2Errors[i]);   // free
		ConvOrders[i] = calloc(NML*NP , sizeof *ConvOrders[i]); // free
	}
	h = calloc(NML*NP , sizeof *h); // free

	// Read in data and compute convergence orders
	for (ML = MLMin; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);

		strcpy(f_name,"../cases/results/");
		strcat(f_name,TestCase); strcat(f_name,"/");
		strcat(f_name,MeshType); strcat(f_name,"/");
		strcat(f_name,"L2errors_");
		sprintf(string,"%dD_",d);   strcat(f_name,string);
		                            strcat(f_name,MeshType);
		sprintf(string,"_ML%d",ML); strcat(f_name,string);
		sprintf(string,"P%d",P);    strcat(f_name,string);
		strcat(f_name,".txt");

		if ((fID = fopen(f_name,"r")) == NULL)
			printf("Error: File: %s, did not open.\n",f_name), exit(1);

		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			i = 0;
			data = StringRead;
			if (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				data += offset;
				h[Indh] = 1.0/pow(tmp_d,1.0/d);
			}
			while (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				L2Errors[i++][Indh] = tmp_d;
				data += offset;
			}
		}
		fclose(fID);
	}}

	for (ML = MLMin+1; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);
		for (i = 0; i < NVars; i++) {
			if (fabs(L2Errors[i][Indh]) > EPS && fabs(L2Errors[i][Indh-NP]) > EPS)
				ConvOrders[i][Indh] = log10(L2Errors[i][Indh]/L2Errors[i][Indh-NP])/log10(h[Indh]/h[Indh-NP]);
		}
	}}

	*pass = 1;

	ML = MLMax;
	for (i = 0; i < NVars; i++) {
		if (!VarsToCheck[i])
			continue;

		for (P = PMin; P <= PMax; P++) {
			Indh = (ML-MLMin)*NP+(P-PMin);
			if ((ConvOrders[i][Indh]-(P+OrderIncrement[i])) < -0.125) {
				*pass = 0;

				printf("i = %d, P = %d, ConvOrder = (% .3e, % .3e)\n",i,P,ConvOrders[i][Indh],1.0*(P+OrderIncrement[i]));
				break;
			}
		}
	}

	free(VarsToCheck);
	free(OrderIncrement);

	if (!(*pass)) {
printf("Re-enable printing here.\n");
//		for (i = 0; i < NVars; i++)
//			array_print_d(NML,NP,ConvOrders[i],'R');
	} else {
		TestDB.Npass++;
	}
printf("ViscousFluxType: %d\n",DB.ViscousFluxType);
printf("h:\n");
array_print_d(NML,NP,h,'R');
printf("L2Errors: \n");
for (i = 0; i < NVars; i++)
array_print_d(NML,NP,L2Errors[i],'R');
printf("Conv Orders: \n");
for (i = 0; i < NVars; i++)
array_print_d(NML,NP,ConvOrders[i],'R');

	for (i = 0; i < NVars; i++) {
		free(L2Errors[i]);
		free(ConvOrders[i]);
	}
	free(L2Errors);
	free(ConvOrders);
	free(h);
}

static void compute_normal(const double *XYZ1, const double *XYZ2, const double *XYZ3, double *n)
{
	/*
	 *	Purpose:
	 *		Compute normal vector to a plane defined by three points.
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
}

void evaluate_mesh_regularity(void)
{
	/*
	 *	Purpose:
	 *		Compute the ratio of enclosing to enclosed spheres by each TET ELEMENT and return the maximum.
	 *
	 *	Comments:
	 *		The mesh regularity check is only necessary for TET refinement as all other ELEMENTs refine into pieces
	 *		having the same shape but with all edge lengths scaled by 1/2 when isotropic h-refinement is performed.
	 *
	 *		The exact radii of the maximal radius enclosed and minimal radius enclosing spheres are calculated as
	 *		follows:
	 *
	 *		Enclosed:
	 *			The enclosed sphere must have a point on each FACET with a vector of length r in the direction normal
	 *			to the surface which starts at this point and ends at the sphere center. Defining this point based on
	 *			the two independent barycentric coordinates on the FACET and noting the four equations arising from the
	 *			four FACETs, a system of 12 equations (4 equations * 3 coordinates) with 12 unknowns (4 * 2 barycentric
	 *			coordinates + 3 centroid coordinates + 1 radius) can be solved for the radius. Surface normal vectors
	 *			are computed based on the equation of the plane describing each FACET.
	 *
	 *		Enclosing:
	 *			There are three possibilities to consider for this case.
	 *
	 *			1) Two   corners of the TET touch the sphere.
	 *			2) Three corners of the TET touch the sphere.
	 *			3) Four  corners of the TET touch the sphere (circumsphere).
	 *
	 *			The respective radii for the cases above are determined based on the vertices of the longest line, the
	 *			triangle formed from the three longest edges, and the circumsphere.
	 *
	 *			Algorithm:
	 *				Check if all four corners fit within the sphere formed from (in the following order):
	 *				- Sphere  of radius (lmax/2) centered at the midpoint of the longest line;
	 *				- Spheres of radius determined from the triangles formed from each face in asending order of radius;
	 *					In this case, the sphere equation is determined from the three coordinates and from the fact
	 *					that the sphere center must lie in the plane of the three vertices.
	 *
	 *				As soon as a condition is met, the appropriate sphere has been found, if neither of the above
	 *				two conditions are met, the circumsphere provides the correct radius.
	 *
	 */

	unsigned int i, j, k, dim, f, d, e, c, Nf, Ne, Nc, iMax, jMax, *IndsE, *IndsF, flag;
	int          **piv;
	double       r_ratio, r, rIn, rOut, *XYZ, *ones, *n, *nNorm, **LHS, **RHS, *lenE, *XYZc, *abcF, *rF, *XYZcT;

	struct S_VOLUME *VOLUME;

	d  = 3;
	Nf = 4;
	Ne = 6;
	Nc = 4;

	unsigned int IndF[3*4]  = { 1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2 },
	             IndE[2*6]  = { 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3 };

	XYZ = malloc(Nc*d * sizeof *XYZ); // free

	ones = malloc(d * sizeof *ones); // free
	for (dim = 0; dim < d; dim++)
		ones[dim] = 1.0;

	n     = malloc(Nf*d * sizeof *n);     // free
	nNorm = malloc(Nf   * sizeof *nNorm); // free

	LHS = malloc(3 * sizeof *LHS); // free
	RHS = malloc(3 * sizeof *RHS); // free
	piv = malloc(3 * sizeof *piv); // free

	LHS[0] = malloc(d*d         * sizeof *LHS[0]); // free
	LHS[1] = malloc(Nf*d*Nf*d   * sizeof *LHS[1]); // free
	LHS[2] = malloc((d+1)*(d+1) * sizeof *LHS[2]); // free
	RHS[0] = malloc(d*1         * sizeof *RHS[0]); // free
	RHS[1] = malloc(Nf*d*1      * sizeof *RHS[1]); // free
	RHS[2] = malloc((d+1)*1     * sizeof *RHS[2]); // free
	piv[0] = malloc(d*1         * sizeof *piv[0]); // free
	piv[1] = malloc(Nf*d*1      * sizeof *piv[1]); // free
	piv[2] = malloc((d+1)*1     * sizeof *piv[2]); // free

	IndsE = malloc(Ne  * sizeof *IndsE); // free
	IndsF = malloc(Nf  * sizeof *IndsF); // free
	lenE  = malloc(Ne  * sizeof *lenE);  // free
	XYZc  = malloc(d   * sizeof *XYZc);  // free
	XYZcT = malloc(d   * sizeof *XYZcT); // free

	abcF = malloc(Nf*d * sizeof *abcF); // free
	rF   = malloc(Nf*1 * sizeof *rF);   // free

	r_ratio = 1.0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->type != TET)
			continue;

		// Obtain vertex coordinates (Row-major)
		for (i = 0; i < Nc; i++) {
		for (j = 0; j < d; j++) {
			XYZ[i*d+j] = VOLUME->XYZ_vC[j*Nc+i];
		}}
array_print_d(Nc,d,XYZ,'R');

		// XYZ coordinates of TET center
		for (j = 0; j < d; j++) {
			XYZcT[j] = 0.0;
			for (i = 0; i < Nc; i++)
				XYZcT[j] += XYZ[i*d+j];
			XYZcT[j] /= 1.0*Nc;
		}

		// Evaluate radius of enclosed sphere

		// 1) Normal vector computation
		for (f = 0; f < Nf; f++) {
			// Cross-product
			compute_normal(&XYZ[IndF[f*d+0]*d],&XYZ[IndF[f*d+1]*d],&XYZ[IndF[f*d+2]*d],RHS[0]);

			// Ensure that normal points outwards
			for (i = 0; i < d; i++) {
				XYZc[i] = 0.0;
				for (j = 0; j < d; j++)
					XYZc[i] += XYZ[IndF[f*d+j]*d+i];
				XYZc[i] /= 3.0;
			}

			r = array_norm_diff_d(d,XYZc,XYZcT,"L2");

			for (i = 0; i < d; i++)
				XYZc[i] += RHS[0][i]*1e1*EPS;

			if (array_norm_diff_d(d,XYZc,XYZcT,"L2")-r < 0.0) {
				for (i = 0; i < d; i++)
					RHS[0][i] *= -1.0;
			}

			nNorm[f] = array_norm_d(d,RHS[0],"L2");
			for (i = 0; i < d; i++)
				n[f*d+i] = RHS[0][i]/nNorm[f];
		}
array_print_d(Nf,d,n,'R');

		// 2) Assemble LHS, RHS ans solve system to obtain rIn
		for (i = 0, iMax = Nf*d*Nf*d; i < iMax; i++)
			LHS[1][i] = 0.0;

		jMax = Nf*d;
		for (f = 0; f < Nf; f++) {
			for (i = f*d; i < (f+1)*d; i++) {
				for (j = f*(d-1); j < (f+1)*(d-1); j++) {
					LHS[1][i*jMax+j] = XYZ[IndF[f*d+(j%(d-1))]*d+(i%d)]-XYZ[IndF[(f+1)*d-1]*d+(i%d)];
				}
				LHS[1][i*jMax+8]       = -n[f*d+(i%d)];
				LHS[1][i*jMax+9+(i%d)] = -1.0;

				if (f < Nf-1)
					RHS[1][i] = -XYZ[3*d+(i%d)];
				else
					RHS[1][i] = -XYZ[2*d+(i%d)];
			}
		}

		if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,Nf*d,1,LHS[1],Nf*d,piv[1],RHS[1],1) > 0)
			printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

		rIn = RHS[1][8];


		// Evaluate radius of enclosing sphere

		// 1) Compute length of all edges
		for (e = 0; e < Ne; e++) {
			lenE[e]  = array_norm_diff_d(d,&XYZ[IndE[e*2+0]*d],&XYZ[IndE[e*2+1]*d],"L2");
			IndsE[e] = e;
		}
		PetscSortRealWithPermutation(Ne,lenE,(int *) IndsE);

		// Sphere radius for case 1.
		rOut = 0.5*lenE[IndsE[Ne-1]];

		for (k = 0; k < d; k++) {
			if (k == 0) {
				// Case 1: 2 Corners touch the sphere

				// Compute sphere center
				for (j = 0; j < d; j++)
					XYZc[j] = 0.5*(XYZ[IndE[IndsE[Ne-1]*2+0]*d+j]+XYZ[IndE[IndsE[Ne-1]*2+1]*d+j]);

				// Check if all corners are enclosed by the sphere
				flag = 1;
				for (c = 0; c < Nc; c++) {
					r = array_norm_diff_d(d,&XYZ[c*d],XYZc,"L2");
					if (r > rOut+EPS)
						flag = 0;
				}
				printf("%d\n",flag);

//				printf("rOut (%d) % .3e\n",k,rOut);
				if (flag)
					break;
			} else if (k == 1) {
				// Case 2: 3 Corners touch the sphere

				// Find the spheres generated by each face and check in order of ascending radius
				jMax = d+1;
				for (f = 0; f < Nf; f++) {
					// Equation of the plane for each face already determined above (coefs == n*nNorm)

					// Solve for sphere center and radius
					for (i = 0; i < d; i++) {
						for (j = 0; j < d; j++)
							LHS[2][i*jMax+j] = 2.0*XYZ[IndF[f*d+i]*d+j];
						LHS[2][i*jMax+d] = 1.0;
						RHS[2][i] = pow(array_norm_d(d,&XYZ[IndF[f*d+i]*d],"L2"),2.0);
					}
					i = d;
					for (j = 0; j < d; j++)
						LHS[2][i*jMax+j] = n[f*d+j]*nNorm[f];
					LHS[2][i*jMax+d] = 0.0;
// need to adjust this '1.0' here, should be a*x3+b*y3+c*z3
					RHS[2][i] = 1.0;

					if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
						printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

					for (i = 0; i < d; i++)
						abcF[f*d+i] = RHS[2][i];
					rF[f] = sqrt(RHS[2][3]+pow(array_norm_d(d,RHS[2],"L2"),2.0));

					IndsF[f] = f;
				}
				PetscSortRealWithPermutation(Nf,rF,(int *) IndsF);

				for (f = 0; f < Nf; f++) {
					XYZc = &abcF[IndsF[f]*d];
					rOut = rF[IndsF[f]];

					// Check if all corners are enclosed by the sphere
					flag = 1;
					for (c = 0; c < Nc; c++) {
						r = array_norm_diff_d(d,&XYZ[c*d],XYZc,"L2");
						if (r > rOut+EPS)
							flag = 0;
					}

//					printf("rOut (%d, %d) % .3e\n",k,f,rOut);
					if (flag)
						break;
				}
			} else if (k == 2) {
				// Case 3: 4 Corners touch the sphere

				// Compute sphere radius
				for (i = 0; i < d+1; i++) {
					for (j = 0; j < d; j++)
						LHS[2][i*jMax+j] = 2.0*XYZ[i*d+j];
					LHS[2][i*jMax+d] = 1.0;
					RHS[2][i] = pow(array_norm_d(d,&XYZ[i*d],"L2"),2.0);
				}

				if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
					printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

				rOut = sqrt(RHS[2][3]+pow(array_norm_d(d,RHS[2],"L2"),2.0));
//				printf("rOut (%d) % .3e\n",k,rOut);
			}
		}

		if (rOut/rIn > r_ratio)
			r_ratio = rOut/rIn;
printf("%d % .3e % .3e % .3e % .3e\n",k,r_ratio,rOut,rIn,sqrt(0.375)*lenE[IndsE[Ne-1]]);

		if (rOut > sqrt(0.375)*lenE[IndsE[Ne-1]]+EPS)
			printf("Error: rOut larger than upper bound on sphere radius.\n"), EXIT_MSG;
	}
EXIT_MSG;
	free(XYZ);
	free(ones);
	free(n);
	free(nNorm);

	array_free2_d(3,LHS);
	array_free2_d(3,RHS);
	array_free2_i(3,piv);

	free(IndsE);
	free(lenE);
	free(XYZc);
	free(XYZcT);

	free(abcF);
	free(rF);
}
