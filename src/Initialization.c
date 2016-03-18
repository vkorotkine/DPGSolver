#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petscksp.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

// MODIFY THE COMMENTS 

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
      Testing              : Run tests for help debugging.
                             Options : 0, 1

  References:

*/

struct S_DB DB;

void Initialization(int nargc, char *argv[]) {
  // Check for presence of '.ctrl' file name input
  char *TestCase;
  TestCase  = malloc(STRLEN_MAX * sizeof *TestCase); // keep

  if (nargc >= 2) TestCase = argv[1];
  else            { printf("Please supply prefix in the compile command.\n"); exit(1); }

  // Set DB Parameters  
  //DB.t_par      = 0; // ToBeModified (Likely initialize all times needed here)

  char *MeshType, *Form, *NodeType, *BasisType, *MeshFile;
  char   *control_file, *string, *dummyc;
  FILE   *fID;

  control_file = malloc(STRLEN_MAX * sizeof *control_file); // tbd
  string       = malloc(STRLEN_MAX * sizeof *string); // tbd
  dummyc       = malloc(STRLEN_MAX * sizeof *dummyc); // tbd

  MeshType  = malloc(STRLEN_MIN * sizeof *MeshType); // keep
  Form      = malloc(STRLEN_MIN * sizeof *Form); // keep
  NodeType  = malloc(STRLEN_MIN * sizeof *NodeType); // keep
  BasisType = malloc(STRLEN_MIN * sizeof *BasisType); // keep
  MeshFile  = malloc(STRLEN_MAX * sizeof *MeshFile); // keep

  // Open control file
  strcpy(control_file,TestCase);
  strcat(control_file,".ctrl");

  fID = fopen(control_file,"r");

  if (fID == NULL) { printf("Control file: %s not present.\n",control_file); exit(1); }

  fscanf(fID,"%[^\n]\n",&string[0]);

  while(!feof(fID))
  {
    fscanf(fID,"%[^\n]\n",&string[0]);
    
    if (strstr(string,"Dimension")  != NULL) sscanf(string,"%s %d",dummyc,&DB.d);
    if (strstr(string,"ML")         != NULL) sscanf(string,"%s %d",dummyc,&DB.ML);
    if (strstr(string,"MeshType")   != NULL) sscanf(string,"%s %s",dummyc,MeshType);
    if (strstr(string,"Form")       != NULL) sscanf(string,"%s %s",dummyc,Form);
    if (strstr(string,"NodeType")   != NULL) sscanf(string,"%s %s",dummyc,NodeType);
    if (strstr(string,"BasisType")  != NULL) sscanf(string,"%s %s",dummyc,BasisType);
    if (strstr(string,"Vectorized") != NULL) sscanf(string,"%s %d",dummyc,&DB.Vectorized);
    if (strstr(string,"EFE")        != NULL) sscanf(string,"%s %d",dummyc,&DB.EFE);
    if (strstr(string,"Collocated") != NULL) sscanf(string,"%s %d",dummyc,&DB.Collocated);
    if (strstr(string,"Adaptive")   != NULL) sscanf(string,"%s %d",dummyc,&DB.Adaptive);
    if (strstr(string,"PGlobal")    != NULL) sscanf(string,"%s %d",dummyc,&DB.P);
    if (strstr(string,"PMax")       != NULL) sscanf(string,"%s %d",dummyc,&DB.PMax);
    if (strstr(string,"Restart")    != NULL) sscanf(string,"%s %d",dummyc,&DB.Restart);
    if (strstr(string,"Testing")    != NULL) sscanf(string,"%s %d",dummyc,&DB.Testing);

    // Mesh file
    if (strstr(string,"BEGIN MESH") != NULL)
    {
      // Note: Scanning only the path here.
      char MeshPath[STRLEN_MAX];
      fscanf(fID,"%s %s\n",dummyc,MeshPath);

      // Note: sprintf prints appends a null character to the character sequence.
      char d[10], ML[10];
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

  // Initialize DB Parameters
  DB.TestCase  = TestCase;
  DB.MeshType  = MeshType;
  DB.Form      = Form;
  DB.NodeType  = NodeType;
  DB.BasisType = BasisType;
  DB.MeshFile  = MeshFile;

  // Print some information
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
