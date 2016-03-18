#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>

#include "database.h"
//#include "parameters.h"
#include "functions.h"

/* 
  Purpose:
    Allocate memory and initialize new structures.
 
  Comments:
  
  Notation:

  Memory freed:

  References:

*/

struct S_ELEMENT *New_ELEMENT(void) {
  int i, iMax;

  struct S_ELEMENT *ELEMENT;
  ELEMENT = malloc(sizeof *ELEMENT); // free

  ELEMENT->present = 0;
  ELEMENT->type    = -1;
  ELEMENT->d       = -1;
  ELEMENT->Nve     = 0;
  ELEMENT->Nf      = 0;

  ELEMENT->Nfve = malloc(2    * sizeof *(ELEMENT->Nfve)); // free
  ELEMENT->VeC  = malloc(8    * sizeof *(ELEMENT->VeC)); // free
  ELEMENT->VeE  = malloc(12*2 * sizeof *(ELEMENT->VeE)); // free
  ELEMENT->VeF  = malloc(6*4  * sizeof *(ELEMENT->VeF)); // free

  for (i = 0; i < 2; i++)                 ELEMENT->Nfve[i] = -1;
  for (i = 0; i < 8; i++)                 ELEMENT->VeC[i]  = -1;
  for (i = 0, iMax = 12*2; i < iMax; i++) ELEMENT->VeE[i]  = -1;
  for (i = 0, iMax = 6*4;  i < iMax; i++) ELEMENT->VeF[i]  = -1;

  ELEMENT->next = NULL;

  return ELEMENT;
}
