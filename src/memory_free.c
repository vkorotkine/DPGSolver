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

void memory_free(void) {
  // DB Parameters

    // Initialization
    free(DB.TestCase), free(DB.MeshType), free(DB.Form), free(DB.NodeType), free(DB.BasisType), free(DB.MeshFile);

    // Preprocessing

      // SetupParameters
      free(DB.Parametrization);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeS);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeF);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeFrs);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeFrc);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeIfs);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeIfc);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeIvs);
      array_free3_c(DB.NP,DB.NDE,DB.NodeTypeIvc);

      free(DB.PGc);
      free(DB.PF);
      array_free2_i(DB.NP,DB.SF_BE);
      array_free2_i(DB.NP,DB.PCs);
      array_free2_i(DB.NP,DB.PCc);
      array_free2_i(DB.NP,DB.PJs);
      array_free2_i(DB.NP,DB.PJc);
      array_free2_i(DB.NP,DB.PFrs);
      array_free2_i(DB.NP,DB.PFrc);
      array_free2_i(DB.NP,DB.PIfs);
      array_free2_i(DB.NP,DB.PIfc);
      array_free2_i(DB.NP,DB.PIvs);
      array_free2_i(DB.NP,DB.PIvc);

      // SetupMesh
      free(DB.PVe), free(DB.NE), free(DB.EType), free(DB.ETags), free(DB.EToVe), free(DB.EToPrt);
      free(DB.VToV), free(DB.VToF), free(DB.VToGF), free(DB.VToBC), free(DB.GFToVe), free(DB.VC), free(DB.GFC);
      free(DB.VeXYZ);

  // ELEMENTs
  struct S_ELEMENT *ELEMENT, *ELEMENTnext;

  ELEMENT = DB.ELEMENT;
  do {
    ELEMENTnext = ELEMENT->next;
    memory_destructor_E(ELEMENT);
    ELEMENT = ELEMENTnext;
  } while(ELEMENTnext != NULL);


}
