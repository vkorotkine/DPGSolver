// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "initialization.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

/*
 *	Purpose:
 *		Declare global variables and read in parameters from the '.ctrl' file.
 *
 *	Comments:
 *
 *	Notation:
 *
 *		Code Parameters
 *			PDE           : PDE to be solved.
 *			PDESpecifier  : Specifier indicating restrictions on the PDE. Also used to choose appropriate mesh based on
 *			                corresponding boundary conditions.
 *
 *			Geometry      : Geometry on which the case is run.
 *			GeomSpecifier : Specifier indicating restrictions on the Geometry. Used to choose appropriate mesh based on
 *			                variable geometry parameters.
 *
 *			Dimension    : d = 2-3. Note that 1D is not supported for any cases and that functions supporting 1D are
 *			               used either for testing purposes or to construct higher dimensional variants.
 *
 *			MeshPath   : Path to the meshes root directory.
 *			MeshType   : Specifies element types used and if the mesh is 'ToBeCurved' in setup_geometry.c
 *			MeshLevel  : (M)esh (L)evel - Uniform refinement level from base mesh (ML = 0)
 *
 *			Form       : Form of the equations (i.e. how many times they are itnegrated by parts)
 *			             Options: Weak
 *			                      Strong
 *			NodeType   : Type of VOLUME nodes to use for different element types
 *			             Options: (T)ensor(P)roduct : (G)auss(L)egendre
 *			                                          (G)auss(L)obatto(L)egendre
 *			                      (SI)mplex         : (A)lpha(O)ptimized
 *			                                          (EQ)ui-spaced
 *			                                          (W)illiams(S)hun - 2D
 *			                                          (S)hun(H)am      - 3D => Denoted as WSH in general
 *			                                          (W)itherden(V)incent
 *			                      WEDGE             : Combination of TP and SI nodes
 *			                      (PYR)amid         : (G)auss(L)egendre
 *			                                          (G)auss(L)obatto(L)egendre
 *			                                          (G)auss(J)acobi
 *			                                          (W)itherden(V)incent
 *			                                          (W)itherden(V)incent(H)exTo(P)YR
 *			BasisType  : Type of basis functions
 *			             Options: Nodal
 *			                      Modal
 *			Vectorized : Type of vectorization to use (ToBeModified)
 *			             Options: 0 (None)
 *			                      1 (Elements grouped by type/order)
 *			EFE        : (E)exact (F)lux (E)valuation - Reduces aliasing if enabled (The analogue in the strong form is
 *			             the (C)hain(R)ule approach)
 *			             Options: 0 (Not used)
 *			                      1 (Used)
 *			Collocated : Specify whether VOLUME nodes should be collocated (Solution/Flux/Flux in reference
 *			             space/Integration nodes)
 *			             Options: 0 (Not collocated)
 *			                      1 (Collocated)
 *			Adapt      : Specify whether adaptation should be used and, if yes, of which type
 *			             Options: 0 (None)
 *			                      1 (p)
 *			                      2 (h)
 *			                      3 (hp)
 *
 *			PGlobal    : Default global (P)olynomial order to be used.
 *			PMax       : (Max)imum (p)olynomial order to be used.
 *			LevelsMax  : (Max)imum number of h refinement (levels) permitted.
 *
 *			Restart    : Specify whether the solution initialization should be based on a previous solution. Currently
 *			             unsupported. (ToBeModified)
 *			             Options: -1           (None)
 *			                       0           (Restart based on solution of order P-1)
 *			                       Iteration # (Restart based on solution of order P at specified iteration #)
 *
 *			Testing    : Run tests for standard checks.
 *			             Options: 0 (No testing)
 *			                      1 (Testing)
 *
 *			ToBeDeleted (Move to initialize_test_case):
 *			BumpOrder  : Specifies polynomial order of "Bumps" along the lower surface of the domain for the
 *			             PolynomialBump test case. The two orders correspond to the two adjacent regions moving away
 *			             from the origin.
 *			             Options: {2 0}, {2 1}, {2 2}, {4 0}
 *			For the PolynomialBump test case, it is advantageous to use BumpOrder = {2,2} as opposed to {4,0} despite
 *			each of them having identical C1 smoothness because the arc length can be computed analytically for
 *			polynomials of order 2 or less, which is much faster.
 *
 *	References:
 *
 */

void initialization(int nargc, char **argv)
{
	// Set DB Parameters
	//DB.t_par      = 0; // ToBeModified (Likely initialize all times needed here)

	// Check for presence of '.ctrl' file name input
	DB.TestCase = malloc(STRLEN_MAX * sizeof *(DB.TestCase)); // keep

	if (nargc >= 2)
		strcpy(DB.TestCase,argv[1]);
	else
		printf("Error: Prefix is absent in the compile command.\n"), EXIT_MSG;

	DB.PDE           = calloc(STRLEN_MAX , sizeof *(DB.PDE));           // keep
	DB.PDESpecifier  = calloc(STRLEN_MAX , sizeof *(DB.PDESpecifier));  // keep
	DB.Geometry      = calloc(STRLEN_MAX , sizeof *(DB.Geometry));      // keep
	DB.GeomSpecifier = calloc(STRLEN_MAX , sizeof *(DB.GeomSpecifier)); // keep

	DB.MeshPath      = calloc(STRLEN_MAX , sizeof *(DB.MeshPath));      // keep
	DB.MeshType      = calloc(STRLEN_MIN , sizeof *(DB.MeshType));      // keep

	char *MeshCurving;
	MeshCurving      = calloc(STRLEN_MIN , sizeof *(MeshCurving));      // free

	DB.Form      = calloc(STRLEN_MIN , sizeof *(DB.Form));      // keep
	DB.NodeType  = calloc(STRLEN_MIN , sizeof *(DB.NodeType));  // keep
	DB.BasisType = calloc(STRLEN_MIN , sizeof *(DB.BasisType)); // keep
	DB.BumpOrder = calloc(2          , sizeof *(DB.BumpOrder)); // keep


	// Open control file
	char *ControlFile = malloc(STRLEN_MAX * sizeof *ControlFile); // free
	strcpy(ControlFile,"control_files/");
	if (strstr(DB.TestCase,"Euler"))
		strcat(ControlFile,"main/Euler/");
	else if (strstr(DB.TestCase,"Test"))
		; // Do nothing
	else
		EXIT_UNSUPPORTED;

	strcat(ControlFile,DB.TestCase);
	strcat(ControlFile,".ctrl");

	FILE *fID;
	if ((fID = fopen(ControlFile,"r")) == NULL)
		printf("Error: Control file: %s not present.\n",ControlFile), EXIT_MSG;
	free(ControlFile);

	// Read input information
	char *StringRead, *dummys;
	StringRead  = malloc(STRLEN_MAX * sizeof *StringRead);  // free
	dummys      = malloc(STRLEN_MAX * sizeof *dummys);      // free

	while(fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		if (strstr(StringRead,"PDEName"))       sscanf(StringRead,"%s %s",dummys,DB.PDE);
		if (strstr(StringRead,"PDESpecifier"))  sscanf(StringRead,"%s %s",dummys,DB.PDESpecifier);

		if (strstr(StringRead,"Geometry"))      sscanf(StringRead,"%s %s",dummys,DB.Geometry);
		if (strstr(StringRead,"GeomSpecifier")) sscanf(StringRead,"%s %s",dummys,DB.GeomSpecifier);

		if (strstr(StringRead,"Dimension"))     sscanf(StringRead,"%s %d",dummys,&DB.d);

		if (strstr(StringRead,"MeshPath"))      sscanf(StringRead,"%s %s",dummys,DB.MeshPath);
		if (strstr(StringRead,"MeshType"))      sscanf(StringRead,"%s %s",dummys,DB.MeshType);
		if (strstr(StringRead,"MeshCurving"))   sscanf(StringRead,"%s %s",dummys,MeshCurving);
		if (strstr(StringRead,"MeshLevel"))     sscanf(StringRead,"%s %d",dummys,&DB.ML);

		if (strstr(StringRead,"Form"))       sscanf(StringRead,"%s %s",dummys,DB.Form);
		if (strstr(StringRead,"NodeType"))   sscanf(StringRead,"%s %s",dummys,DB.NodeType);
		if (strstr(StringRead,"BasisType"))  sscanf(StringRead,"%s %s",dummys,DB.BasisType);
		if (strstr(StringRead,"Vectorized")) sscanf(StringRead,"%s %d",dummys,&DB.Vectorized);
		if (strstr(StringRead,"EFE"))        sscanf(StringRead,"%s %d",dummys,&DB.EFE);
		if (strstr(StringRead,"Collocated")) sscanf(StringRead,"%s %d",dummys,&DB.Collocated);
		if (strstr(StringRead,"Adapt"))      sscanf(StringRead,"%s %d",dummys,&DB.Adapt);
		if (strstr(StringRead,"PGlobal"))    sscanf(StringRead,"%s %d",dummys,&DB.PGlobal);
		if (strstr(StringRead,"PMax"))       sscanf(StringRead,"%s %d",dummys,&DB.PMax);
		if (strstr(StringRead,"LevelsMax"))  sscanf(StringRead,"%s %d",dummys,&DB.LevelsMax);
		if (strstr(StringRead,"Restart"))    sscanf(StringRead,"%s %d",dummys,&DB.Restart);
		if (strstr(StringRead,"Testing"))    sscanf(StringRead,"%s %d",dummys,&DB.Testing);

		if (strstr(StringRead,"BumpOrder")) {
			EXIT_UNSUPPORTED;
			// The BumpOrder should be specified as part the the GeomSpecifier and used to set the appropriate
			// parameters in initialize_test_cases. ToBeDeleted
			sscanf(StringRead,"%s %d %d",dummys,&DB.BumpOrder[0],&DB.BumpOrder[1]);
		}
	}
	free(StringRead);
	free(dummys);
	fclose(fID);

	if (!strstr(MeshCurving,"Straight")) {
		strcat(MeshCurving,DB.MeshType);
		strcpy(DB.MeshType,MeshCurving);
	}
	free(MeshCurving);

	// Set up MeshFile
	DB.MeshFile = calloc(STRLEN_MAX , sizeof *(DB.MeshFile)); // keep

	char *d, *ML;
	d  = malloc(STRLEN_MIN * sizeof *d);  // free
	ML = malloc(STRLEN_MIN * sizeof *ML); // free

	sprintf(d,"%d",DB.d);
	sprintf(ML,"%d",DB.ML);

	strcpy(DB.MeshFile,"");

	strcat(DB.MeshFile,DB.MeshPath);
	strcat(DB.MeshFile,DB.Geometry); strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.PDE);      strcat(DB.MeshFile,"/");
	if (DB.PDESpecifier)
		strcat(DB.MeshFile,DB.PDESpecifier);  strcat(DB.MeshFile,"/");
	if (DB.GeomSpecifier)
		strcat(DB.MeshFile,DB.GeomSpecifier); strcat(DB.MeshFile,"/");
	strcat(DB.MeshFile,DB.Geometry);
	strcat(DB.MeshFile,strcat(d,"D_"));
	strcat(DB.MeshFile,DB.MeshType);
	strcat(DB.MeshFile,strcat(ML,"x.msh"));

	free(d);
	free(ML);

	// Print some information
	if (!DB.MPIrank && !DB.Testing) {
		printf("\n\nRunning the %s test case using the %s mesh type in %dD on mesh level %d.\n\n",
		       DB.TestCase,DB.MeshType,DB.d,DB.ML);
		printf("Parameters:\n\n");
		printf("Form       : %s\n",    DB.Form);
		printf("NodeType   : %s\n",    DB.NodeType);
		printf("BasisType  : %s\n\n",  DB.BasisType);
		printf("Vectorized : %d\n",    DB.Vectorized);
		printf("EFE        : %d\n",    DB.EFE);
		printf("Collocated : %d\n",    DB.Collocated);
		printf("Adapt      : %d\n\n",  DB.Adapt);
		printf("P          : %d\n",    DB.PGlobal);
		printf("PMax       : %d\n",    DB.PMax);
		printf("LevelsMax  : %d\n\n\n",DB.LevelsMax);
	}

	if ((DB.Adapt == ADAPT_H || DB.Adapt == ADAPT_HP) && DB.LevelsMax == 0)
		printf("\n\n***** Warning: LevelsMax should be greater than 0 if h-adaptation is enabled. *****\n\n\n");
}
