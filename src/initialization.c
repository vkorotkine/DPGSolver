// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"

/*
 *	Purpose:
 *		Declare global variables and read in parameters from the '.ctrl' file.
 *
 *	Comments:
 *		For the PolynomialBump test case, it is advantageous to use BumpOrder = {2,2} as opposed to {4,0} despite each
 *		of them having identical C1 smoothness because the arc length can be computed analytically for polynomials of
 *		order 2 or less, which is much faster.
 *
 *	Notation:
 *
 *		Code Parameters
 *			Dimension  : d = 1-3
 *			ML         : (M)esh (L)evel - Uniform refinement level from base mesh (ML = 0)
 *			MeshType   : Specifies element types used and if the mesh is 'ToBeCurved' in setup_geometry.c
 *			BumpOrder  : Specifies polynomial order of "Bumps" along the lower surface of the domain for the
 *			             PolynomialBump test case. The two orders correspond to the two adjacent regions moving away
 *			             from the origin.
 *			             Options: {2 0}, {2 1}, {2 2}, {4 0}
 *
 *			Form       : Form of the equations (i.e. how many times they are itnegrated by parts)
 *			             Options: Weak
 *			                      Strong
 *			NodeType   : Type of VOLUME nodes to use for different element types
 *			             Options: (T)ensor(P)roduct : (G)auss(L)egendre
 *			                                          (G)auss(L)obatto(L)egendre
 *			                      (SI)mplex         : (A)lpha(O)ptimized
 *			                                          (W)illiams(S)hun - 2D
 *			                                          (S)hun(H)am      - 3D
 *			                      WEDGE             : Combination of TP and SI nodes
 *			                      (PYR)amid         : ToBeModified
 *			BasisType  : Type of basis functions
 *			             Options: Nodal
 *			                      Modal
 *			Vectorized : Type of vectorization to use (ToBeModified)
 *			             Options: 0 (None)
 *			                      1 (Elements grouped by type/order)
 *			EFE        : (E)exact (F)lux (E)valuation - Reduces aliasing if enabled (The analogue in the strong form is
 *			             the CR (Chain-Rule) approach)
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
 *			P          : Polynomial order to be used (not used if p-adaptation is enabled)
 *			PMax       : Maximum polynomial order to be used (used onnly if p-adaptation is enabled)
 *
 *			Restart    : Specify whether the solution initialization should be based on a previous solution
 *			             Options: -1           (None)
 *			                       0           (Restart based on solution of order P-1)
 *			                       Iteration # (Restart based on solution of order P at specified iteration #)
 *
 *			Testing    : Run tests for standard checks.
 *			             Options: 0 (No testing)
 *			                      1 (Testing)
 *
 *	References:
 *
 */

void initialization(int nargc, char **argv)
{
	// Set DB Parameters
	//DB.t_par      = 0; // ToBeModified (Likely initialize all times needed here)

	char         *TestCase, *MeshType, *Form, *NodeType, *BasisType, *MeshFile, *ControlFile, *StringRead, *dummys,
	             *MeshPath, *d, *ML;
	unsigned int *BumpOrder;
	FILE         *fID;

	// Check for presence of '.ctrl' file name input
	TestCase  = malloc(STRLEN_MAX * sizeof *TestCase); // keep

	if (nargc >= 2)
		strcpy(TestCase,argv[1]);
	else
		printf("Error: Prefix is absent in the compile command.\n"), exit(1);

	ControlFile = malloc(STRLEN_MAX * sizeof *ControlFile); // free
	StringRead  = malloc(STRLEN_MAX * sizeof *StringRead);  // free
	dummys      = malloc(STRLEN_MAX * sizeof *dummys);      // free

	MeshType  = malloc(STRLEN_MIN * sizeof *MeshType);  // keep
	Form      = malloc(STRLEN_MIN * sizeof *Form);      // keep
	NodeType  = malloc(STRLEN_MIN * sizeof *NodeType);  // keep
	BasisType = malloc(STRLEN_MIN * sizeof *BasisType); // keep
	MeshFile  = malloc(STRLEN_MAX * sizeof *MeshFile);  // keep

	MeshPath  = malloc(STRLEN_MAX * sizeof *MeshPath); // free
	d         = malloc(STRLEN_MIN * sizeof *d);        // free
	ML        = malloc(STRLEN_MIN * sizeof *ML);       // free

	BumpOrder = calloc(2          , sizeof *BumpOrder); // keep

	// Open control file
	strcpy(ControlFile,TestCase);
	strcat(ControlFile,".ctrl");

	if ((fID = fopen(ControlFile,"r")) == NULL)
		printf("Error: Control file: %s not present.\n",ControlFile), exit(1);
	free(ControlFile);

	while(fscanf(fID,"%[^\n]\n",StringRead) == 1) {

		if (strstr(StringRead,"Dimension")  != NULL) sscanf(StringRead,"%s %d",dummys,&DB.d);
		if (strstr(StringRead,"ML")         != NULL) sscanf(StringRead,"%s %d",dummys,&DB.ML);
		if (strstr(StringRead,"MeshType")   != NULL) sscanf(StringRead,"%s %s",dummys,MeshType);
		if (strstr(StringRead,"BumpOrder")  != NULL) sscanf(StringRead,"%s %d %d",dummys,&BumpOrder[0],&BumpOrder[1]);
		if (strstr(StringRead,"Form")       != NULL) sscanf(StringRead,"%s %s",dummys,Form);
		if (strstr(StringRead,"NodeType")   != NULL) sscanf(StringRead,"%s %s",dummys,NodeType);
		if (strstr(StringRead,"BasisType")  != NULL) sscanf(StringRead,"%s %s",dummys,BasisType);
		if (strstr(StringRead,"Vectorized") != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Vectorized);
		if (strstr(StringRead,"EFE")        != NULL) sscanf(StringRead,"%s %d",dummys,&DB.EFE);
		if (strstr(StringRead,"Collocated") != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Collocated);
		if (strstr(StringRead,"Adapt")      != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Adapt);
		if (strstr(StringRead,"PGlobal")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.PGlobal);
		if (strstr(StringRead,"PMax")       != NULL) sscanf(StringRead,"%s %d",dummys,&DB.PMax);
		if (strstr(StringRead,"Restart")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Restart);
		if (strstr(StringRead,"Testing")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Testing);

		// Mesh file
		if (strstr(StringRead,"BEGIN MESH") != NULL) {
			if(fscanf(fID,"%s %s\n",dummys,MeshPath) == 2) {
				sprintf(d,"%d",DB.d);
				sprintf(ML,"%d",DB.ML);

				strcpy(MeshFile,"");

				strcat(MeshFile,MeshPath);
				strcat(MeshFile,TestCase);
				strcat(MeshFile,"/");
				strcat(MeshFile,TestCase);
				strcat(MeshFile,strcat(d,"D_"));
				strcat(MeshFile,MeshType);
				strcat(MeshFile,strcat(ML,"x.msh"));
			}
		}
	}
	fclose(fID);

	free(StringRead);
	free(dummys);
	free(MeshPath);
	free(d);
	free(ML);

	// Assign DB Parameters
	DB.TestCase  = TestCase;
	DB.MeshType  = MeshType;
	DB.Form      = Form;
	DB.NodeType  = NodeType;
	DB.BasisType = BasisType;
	DB.MeshFile  = MeshFile;
	DB.BumpOrder = BumpOrder;

	// Print some information
	if (!DB.MPIrank) {
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
		printf("Testing    : %d\n\n\n",DB.Testing);
	}
}
