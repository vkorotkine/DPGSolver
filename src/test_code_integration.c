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
#include "setup_Curved.h"
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
	strcat(DB.MeshFile,DB.Geometry);
	strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.Geometry);
	strcat(DB.MeshFile,strcat(d,"D_"));
	strcat(DB.MeshFile,DB.MeshType);
	strcat(DB.MeshFile,strcat(ML,"x.msh"));

	free(d);
	free(ML);
}

static void update_TestCase(void)
{
	if (strstr(DB.TestCase,"Poisson_Ringleb")) {
		strcpy(DB.TestCase,"Poisson");
		strcpy(DB.Geometry,"Ringleb");
	} else if (strstr(DB.TestCase,"Poisson_dm1-Spherical_Section")) {
		strcpy(DB.TestCase,"Poisson");
		strcpy(DB.Geometry,"dm1-Spherical_Section");
	} else if (strstr(DB.TestCase,"L2_proj") ||
	           strstr(DB.TestCase,"update_h")) {
		strcpy(DB.TestCase,"PeriodicVortex");
		strcpy(DB.Geometry,"PeriodicVortex"); // ToBeModified: Rename this.
	} else if (strstr(DB.TestCase,"linearization")) {
		strcpy(DB.TestCase,"SupersonicVortex");
		strcpy(DB.Geometry,"SupersonicVortex"); // ToBeModified: Rename this.
	} else {
		printf("%s\n",DB.TestCase);
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
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
		DB.PGlobal = TestDB.PGlobal;
		if (update_argv == 1)
			DB.ML = TestDB.ML;
		update_TestCase();
		update_MeshFile();
	}

	setup_parameters();
	if (update_argv)
		setup_parameters_L2proj();

	initialize_test_case_parameters();
	setup_mesh();
	setup_operators();
	setup_structures();
	setup_geometry();

	initialize_test_case(Nref);

	if (update_argv == 2)
		mesh_to_level(TestDB.ML);
}

void code_cleanup(void)
{
	mesh_to_level(0);
	memory_free();
}

void code_startup_mod_prmtrs(int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv,
                             const unsigned int phase)
{
	int  MPIrank, MPIsize;

	if (phase == 1) {
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
			DB.PGlobal = TestDB.PGlobal;
			if (update_argv == 1)
				DB.ML      = TestDB.ML;
			update_TestCase();
			update_MeshFile();
		}

		initialize_test_case_parameters();
		setup_parameters();
		if (update_argv)
			setup_parameters_L2proj();
	} else if (phase == 2) {
		setup_mesh();
		setup_operators();
		setup_structures();
		setup_geometry();

		initialize_test_case(Nref);

		if (update_argv == 2)
			mesh_to_level(TestDB.ML);
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}


void check_convergence_orders(const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                              const unsigned int PMax, unsigned int *pass)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
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
printf("ViscousFlux, Blending, Parametrization: %d, %d, %d\n",DB.ViscousFluxType,DB.Blending,DB.Parametrization);
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

void evaluate_mesh_regularity(double *mesh_quality)
{
	/*
	 *	Purpose:
	 *		Compute the ratio of enclosing to enclosed spheres by each TET ELEMENT and return the maximum.
	 *
	 *	Comments:
	 *		The mesh regularity check is only necessary for TET refinement as all other ELEMENTs refine into pieces
	 *		having the same shape but with all edge lengths scaled by 1/2 when isotropic h-refinement is performed.
	 *
	 *		While testing this function, several meshes generated by gmsh returned elements having extremely small
	 *		internal radius (nearly coplanar TET vertices). These elements, being of the worst quality, resulted in a
	 *		growing ratio of the internal to external radii with with the mesh refinement despite obtaining optimal
	 *		convergence orders (L2 norm). As the L2 norm is global, it was thus decided to take the average of ratios
	 *		from each element. Correct results were then obtained (Suboptimal using refine by split, optimal using
	 *		non-nested unstructured meshes).
	 *			-> Potentially look into why gmsh gives such poor quality elements. (ToBeModified)
	 *
	 *		The exact radii of the maximal radius enclosed and minimal radius enclosing spheres are calculated as
	 *		follows:
	 *
	 *		Enclosed:
	 *			Define the planes of each face by a*x+b*y+c*z = d:
	 *				1) Compute [a b c] as the normal vector to the plane and find d by substituting one coordinate.
	 *				2) Normalize [a b c] (and d) and invert the normal if not pointing outwards:
	 *					Given distance from point to plane = |a*x+b*y+c*z-d|/norm([a b c],2), ensure that when
	 *					substituting the coordinate of the vertex not in the plane that the result is negative.
	 *				3) Using the distance formula again:
	 *					r = n1*x_c+n2*y_c+n3*z_c-d (4 eqns and 4 unknowns)
	 *					Solve: [ones(d+1,1) nr]*[r x_c y_c z_c]' = d
	 *
	 *		Enclosing:
	 *			It may be sufficient for the mesh quality measure to simply use the circumsphere radius. (ToBeModified)
	 *
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

	unsigned int i, j, k, dim, f, d, e, c, Nf, Ne, Nc, iMax, jMax, *IndsE, *IndsF, Found, TETcount, NormType;
	int          **piv;
	double       r_ratio, r, rIn, rOut,
	             *XYZ, *XYZdiff, *ones, *n, *nNorm, **LHS, **RHS, *lenE, *XYZc, *abcF, *rF, *XYZcT, *d_p;

	struct S_VOLUME *VOLUME;

	d  = 3;
	Nf = 4;
	Ne = 6;
	Nc = 4;

	unsigned int IndF[3*4]  = { 1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2 },
	             IndE[2*6]  = { 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3 };

	NormType = 0; // Options: 0 (Inf), 2 (L2)

	XYZ     = malloc(Nc*d * sizeof *XYZ);     // free
	XYZdiff = malloc(d    * sizeof *XYZdiff); // free

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
	d_p   = malloc(Nf  * sizeof *d_p);   // free

	abcF = malloc(Nf*d * sizeof *abcF); // free
	rF   = malloc(Nf*1 * sizeof *rF);   // free

	r_ratio  = 0.0;
	TETcount = 0;
	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (VOLUME->type != TET)
			continue;

		// Obtain vertex coordinates (Row-major)
		for (i = 0; i < Nc; i++) {
		for (j = 0; j < d; j++) {
			XYZ[i*d+j] = VOLUME->XYZ_vV[j*Nc+i];
		}}

		// XYZ coordinates of TET center
		for (j = 0; j < d; j++) {
			XYZcT[j] = 0.0;
			for (i = 0; i < Nc; i++)
				XYZcT[j] += XYZ[i*d+j];
			XYZcT[j] /= 1.0*Nc;
		}

		// Evaluate radius of enclosed sphere

		// 1) Normal vector (normalized) computation
		for (f = 0; f < Nf; f++) {
			// Cross-product
			compute_plane(&XYZ[IndF[f*d+0]*d],&XYZ[IndF[f*d+1]*d],&XYZ[IndF[f*d+2]*d],RHS[0],&d_p[f]);

			// Normalize
			nNorm[f] = array_norm_d(d,RHS[0],"L2");
			for (i = 0; i < d; i++)
				n[f*d+i] = RHS[0][i]/nNorm[f];
			d_p[f] /= nNorm[f];

			// Ensure that normal points outwards
			r = -d_p[f];
			for (i = 0; i < d; i++)
				r += n[f*d+i]*XYZ[f*d+i];

			if (r > 0.0) {
				for (i = 0; i < d; i++)
					n[f*d+i] *= -1.0;
				d_p[f] *= -1.0;
			}
		}

		// 2) Assemble LHS, RHS and solve system to obtain rIn
		for (i = 0, iMax = d+1; i < iMax; i++) {
			for (j = 0, jMax = d+1; j < jMax; j++) {
				if (j == 0)
					LHS[2][i*jMax+j] = 1.0;
				else
					LHS[2][i*jMax+j] = n[i*d+(j-1)];
			}
			RHS[2][i] = d_p[i];
		}

		if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
			printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

		rIn = RHS[2][0];


		// Evaluate radius of enclosing sphere

		// 1) Compute length of all edges
		for (e = 0; e < Ne; e++) {
			for (i = 0; i < d; i++)
				XYZdiff[i] = XYZ[IndE[e*2+0]*d+i]-XYZ[IndE[e*2+1]*d+i];
			lenE[e]  = array_norm_d(d,XYZdiff,"L2");
			IndsE[e] = e;
		}
		PetscSortRealWithPermutation(Ne,lenE,(int *) IndsE);

		// Sphere radius for case 1.
		rOut = 0.5*lenE[IndsE[Ne-1]];

		Found = 0;
		for (k = 0; k < d && !Found; k++) {
			if (k == 0) {
				// Case 1: 2 Corners touch the sphere

				// Compute sphere center
				for (j = 0; j < d; j++)
					XYZc[j] = 0.5*(XYZ[IndE[IndsE[Ne-1]*2+0]*d+j]+XYZ[IndE[IndsE[Ne-1]*2+1]*d+j]);

				// Check if all corners are enclosed by the sphere
				Found = 1;
				for (c = 0; c < Nc; c++) {
					for (i = 0; i < d; i++)
						XYZdiff[i] = XYZ[c*d+i]-XYZc[i];
					r = array_norm_d(d,XYZdiff,"L2");
					if (r > rOut+EPS)
						Found = 0;
				}
//				printf("rOut (%d) % .3e\n",k,rOut);
				if (Found)
					break;
			} else if (k == 1) {
				// Case 2: 3 Corners touch the sphere

				// Find the spheres generated by each face and check in order of ascending radius
				jMax = d+1;
				for (f = 0; f < Nf; f++) {
					// Equation of the plane for each face already determined above (coefs == n*nNorm, d == d_p)

					// Solve for sphere center and radius
					for (i = 0; i < d; i++) {
						for (j = 0; j < d; j++)
							LHS[2][i*jMax+j] = 2.0*XYZ[IndF[f*d+i]*d+j];
						LHS[2][i*jMax+d] = 1.0;
						RHS[2][i] = pow(array_norm_d(d,&XYZ[IndF[f*d+i]*d],"L2"),2.0);
					}
					i = d;
					for (j = 0; j < d; j++)
						LHS[2][i*jMax+j] = n[f*d+j];
					LHS[2][i*jMax+d] = 0.0;
					RHS[2][i] = d_p[f];

					if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
						printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

					for (i = 0; i < d; i++)
						abcF[f*d+i] = RHS[2][i];
					rF[f] = sqrt(RHS[2][3]+pow(array_norm_d(d,RHS[2],"L2"),2.0));

					IndsF[f] = f;
				}
				PetscSortRealWithPermutation(Nf,rF,(int *) IndsF);

				for (f = 0; f < Nf; f++) {
					for (i = 0; i < d; i++)
						XYZc[i] = abcF[IndsF[f]*d+i];
					rOut = rF[IndsF[f]];

					// Check if all corners are enclosed by the sphere
					Found = 1;
					for (c = 0; c < Nc; c++) {
						for (i = 0; i < d; i++)
							XYZdiff[i] = XYZ[c*d+i]-XYZc[i];
						r = array_norm_d(d,XYZdiff,"L2");
						if (r-rOut > 1e2*EPS)
							Found = 0;
					}
//					printf("rOut (%d, %d) % .3e\n",k,f,rOut);
					if (Found)
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

		if (NormType == 0) {
			if (rOut/rIn > r_ratio) {
				r_ratio = rOut/rIn;
			}
		} else if (NormType == 2) {
			r_ratio += rOut/rIn;
		}
		TETcount += 1;

//printf("%d % .3e % .3e % .3e % .3e % .3e\n",
//       k,rOut/rIn,rOut,rIn,sqrt(0.375)*lenE[IndsE[Ne-1]],rOut-sqrt(0.375)*lenE[IndsE[Ne-1]]);

		if (rOut - sqrt(0.375)*lenE[IndsE[Ne-1]] > 1e2*EPS)
			printf("Error: rOut larger than upper bound on sphere radius.\n"), EXIT_MSG;
	}
	if (TETcount) {
		if (NormType == 0)
			*mesh_quality = r_ratio;
		else if (NormType == 2)
			*mesh_quality = r_ratio/TETcount;
		else
			printf("Error: Unsupported.\n"), EXIT_MSG;
	} else {
		*mesh_quality = 3.0; // ratio for regular TET
	}

//	printf("%d % .3e\n",DB.ML,*mesh_quality);

	free(XYZ);
	free(XYZdiff);
	free(ones);
	free(n);
	free(nNorm);

	array_free2_d(3,LHS);
	array_free2_d(3,RHS);
	array_free2_i(3,piv);

	free(IndsE);
	free(IndsF);
	free(lenE);
	free(XYZc);
	free(XYZcT);
	free(d_p);

	free(abcF);
	free(rF);
}

void check_mesh_regularity(const double *mesh_quality, const unsigned int NML, unsigned int *pass)
{
	/*
	 *	Purpose:
	 *		Checks that the sequence of refined meshes used satisfy the mesh regularity condition required for optimal
	 *		convergence to be obtained.
	 *
	 *	Comments:
	 *		For the moment, the slopes of the last two intervals are used to assess this.
	 */

	unsigned int i;
	double       slope_quality[2];

	if (NML > 2) {
		for (i = 0; i < 2; i++)
			slope_quality[i] = mesh_quality[NML-2+i]-mesh_quality[NML-3+i];

		if (slope_quality[1] > 0.0) {
			if (slope_quality[1]-slope_quality[0] > 1e-3)
				*pass = 0;
		}
		printf("\nMesh quality:");
		array_print_d(1,NML,mesh_quality,'R');
	}
}
