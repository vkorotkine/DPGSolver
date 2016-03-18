#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
//#include "functions.h"

/* 
  Purpose:
    Declare global variables and read in parameters from the '.ctrl' file.
 
  Comments:
  
  Notation:
    DB      : D(ata)B(ase)

    Code Parameters
      Dimension            : d = 1-3
      (M)esh(L)evel        : Uniform refinement level from base mesh (ML = 0)
      MeshType             : Specifies element types used and if the mesh is 'ToBeCurved' in SetupGeometry.c

      Form                 : Form of the equations (i.e. how many times they are integrated by parts)
                             Options : Weak, Strong
      NodeType             : Type of VOLUME nodes to use for different element types
                             Options : (T)ensor(P)roduct : (G)auss(L)egendre, (G)auss(L)obatto(L)egendre
                                       (SI)mplex         : (A)lpha(O)ptimized, (W)illiams(S)hun - 2D, (S)hun(H)am - 3D
                                       WEDGE             : Combination of TP and SI nodes
                                       (PYR)amid         : ToBeModified
      BasisType            : Type of basis functions
                             Options : Nodal, Modal
      Vectorized           : Type of vectorization to use
                             Options : 0 (None), 1 (Standard C loops), 2 (BLAS, Dense operators only), 3 (BLAS, All operators)
      (E)exact(F)lux(E)val : Reduce aliasing by using exact flux evaluations (The analogue in the strong form is the
                             (C)hain(R)ule approach)
                             Options : 0, 1
      Collocated           : Specify whether VOLUME nodes should be collocated (Solution/Flux/Flux in reference
                             space/Integration nodes)
                             Options : 0, 1
      Adaptive             : Specify whether adaptation should be used and, if yes, of which type
                             Options : 0 (None), 1 (p), 2 (h), 3 (hp)

      P                    : Polynomial order to be used.         Not used if p-adaptation is enabled
      PMax                 : Maximum polynomial order to be used. Not used unless p-adaptation is enabled

      Restart              : Specify whether the solution initialization should be based on a previous solution
                             Options: -1           (None)
                                       0           (Restart based on solution of order P-1)
                                       Iteration # (Restart based on solution of order P at specified iteration)
      Testing              : Run tests for standard checks.
                             Options : 0, 1

  References:

*/

void Initialization(int nargc, char **argv) {
  // Set DB Parameters  
  //DB.t_par      = 0; // ToBeModified (Likely initialize all times needed here)

  char *TestCase, *MeshType, *Form, *NodeType, *BasisType, *MeshFile, 
       *ControlFile, *StringRead, *dummys, *MeshPath, *d, *ML;
  FILE *fID;

  // Check for presence of '.ctrl' file name input
  TestCase  = malloc(STRLEN_MAX * sizeof *TestCase); // keep

  if (nargc >= 2) strcpy(TestCase,argv[1]);
  else            printf("Please supply prefix in the compile command.\n"), exit(1);

  ControlFile = malloc(STRLEN_MAX * sizeof *ControlFile); // free
  StringRead  = malloc(STRLEN_MAX * sizeof *StringRead); // free
  dummys      = malloc(STRLEN_MAX * sizeof *dummys); // free

  MeshType  = malloc(STRLEN_MIN * sizeof *MeshType); // keep
  Form      = malloc(STRLEN_MIN * sizeof *Form); // keep
  NodeType  = malloc(STRLEN_MIN * sizeof *NodeType); // keep
  BasisType = malloc(STRLEN_MIN * sizeof *BasisType); // keep
  MeshFile  = malloc(STRLEN_MAX * sizeof *MeshFile); // keep

  MeshPath  = malloc(STRLEN_MAX * sizeof *MeshPath); // free
  d         = malloc(STRLEN_MIN * sizeof *d); // free
  ML        = malloc(STRLEN_MIN * sizeof *ML); // free

  // Open control file
  strcpy(ControlFile,TestCase);
  strcat(ControlFile,".ctrl");

  if ((fID = fopen(ControlFile,"r")) == NULL)
    printf("Error: Control file: %s not present.\n",ControlFile), exit(1);
  free(ControlFile);

  fscanf(fID,"%[^\n]\n",StringRead);

  while(!feof(fID)) {
    fscanf(fID,"%[^\n]\n",StringRead);
    
    if (strstr(StringRead,"Dimension")  != NULL) sscanf(StringRead,"%s %d",dummys,&DB.d);
    if (strstr(StringRead,"ML")         != NULL) sscanf(StringRead,"%s %d",dummys,&DB.ML);
    if (strstr(StringRead,"MeshType")   != NULL) sscanf(StringRead,"%s %s",dummys,MeshType);
    if (strstr(StringRead,"Form")       != NULL) sscanf(StringRead,"%s %s",dummys,Form);
    if (strstr(StringRead,"NodeType")   != NULL) sscanf(StringRead,"%s %s",dummys,NodeType);
    if (strstr(StringRead,"BasisType")  != NULL) sscanf(StringRead,"%s %s",dummys,BasisType);
    if (strstr(StringRead,"Vectorized") != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Vectorized);
    if (strstr(StringRead,"EFE")        != NULL) sscanf(StringRead,"%s %d",dummys,&DB.EFE);
    if (strstr(StringRead,"Collocated") != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Collocated);
    if (strstr(StringRead,"Adaptive")   != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Adaptive);
    if (strstr(StringRead,"PGlobal")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.P);
    if (strstr(StringRead,"PMax")       != NULL) sscanf(StringRead,"%s %d",dummys,&DB.PMax);
    if (strstr(StringRead,"Restart")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Restart);
    if (strstr(StringRead,"Testing")    != NULL) sscanf(StringRead,"%s %d",dummys,&DB.Testing);

    // Mesh file
    if (strstr(StringRead,"BEGIN MESH") != NULL) {
      fscanf(fID,"%s %s\n",dummys,MeshPath);

      sprintf(d,"%d",DB.d);
      sprintf(ML,"%d",DB.ML);

      strcpy(MeshFile,"");

      strcat(MeshFile,MeshPath);
      strcat(MeshFile,TestCase); // Add final path argument 
      strcat(MeshFile,"/"); 
      strcat(MeshFile,TestCase);
      strcat(MeshFile,strcat(d,"D_"));
      strcat(MeshFile,MeshType);
      strcat(MeshFile,strcat(ML,"x.msh"));
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

  // Print some information
  if (!DB.MPIrank) {
    printf("\n\nRunning the %s testcase using the %s mesh type in %dD on mesh level %d.\n\n",
           DB.TestCase,DB.MeshType,DB.d,DB.ML);
    printf("Parameters:\n\n");
    printf("Form       : %s\n",    DB.Form);
    printf("NodeType   : %s\n",    DB.NodeType);
    printf("BasisType  : %s\n\n",  DB.BasisType);
    printf("Vectorized : %d\n",    DB.Vectorized);
    printf("EFE        : %d\n",    DB.EFE);
    printf("Collocated : %d\n",    DB.Collocated);
    printf("Adaptive   : %d\n\n",  DB.Adaptive);
    printf("P          : %d\n",    DB.P);
    printf("PMax       : %d\n",    DB.PMax);
    printf("Testing    : %d\n\n\n",DB.Testing);
  }
}
