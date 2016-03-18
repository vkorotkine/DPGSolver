#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

#include "database.h"
//#include "parameters.h"
#include "functions.h"

/* 
  Purpose:
    Free remaining memory.
 
  Comments:
  
  Notation:

  Memory freed:

  References:

*/

struct S_DB DB;

void FreeMemory(void) {

  // DB Parameters

    // Initialization
    free(DB.TestCase), free(DB.MeshType), free(DB.Form), free(DB.NodeType), free(DB.BasisType), free(DB.MeshFile);

    // Preprocessing

      // SetupParameters
      free(DB.Parametrization);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeS);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeF);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeFrs);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeFrc);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeIfs);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeIfc);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeIvs);
      ArrayFree3c(DB.NP,DB.NDE,DB.NodeTypeIvc);

      free(DB.PGc);
      free(DB.PF);
      ArrayFree2i(DB.NP,DB.SF_BE);
      ArrayFree2i(DB.NP,DB.PCs);
      ArrayFree2i(DB.NP,DB.PCc);
      ArrayFree2i(DB.NP,DB.PJs);
      ArrayFree2i(DB.NP,DB.PJc);
      ArrayFree2i(DB.NP,DB.PFrs);
      ArrayFree2i(DB.NP,DB.PFrc);
      ArrayFree2i(DB.NP,DB.PIfs);
      ArrayFree2i(DB.NP,DB.PIfc);
      ArrayFree2i(DB.NP,DB.PIvs);
      ArrayFree2i(DB.NP,DB.PIvc);



}
