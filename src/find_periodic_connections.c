#include <stdlib.h>
#include <stdio.h>

#include "functions.h"

#include "petscsys.h"

/*
*	Purpose:
*	Given a list of periodic vertex correspondence, find all possible matches between the vertices.
*
*	Comments:
*	Each row of the list is sorted in ascending order. Hence, reversed entries are redundant and not included.
*
*	Notation:
*	PVe : List of periodic vertices (pve x 2 array)
*	pve : Index of last row in PVe
*
*	Example:
*
*	Input:
*	PVe = 
*	1 0
*	2 8
*	1 3
*	8 3
*	3 2
*	8 1
*	5 8
**
*	pve = 7
*
*	Pass 1:
*
*	PVeMatches:
*	1 0 0 0 0 0 0 0
*	0 3 8 0 0 0 0 0
*	3 8 0 0 0 0 0 0
**	1 2 8 0 0 0 0 0
*	8 0 0 0 0 0 0 0
*	1 2 3 5 0 0 0 0
*
*	Updated PVe (Sorted):
*	0 1
*	0 3
*	0 8
*	1 2
*	1 3
*	1 5
**	1 8
*	2 3
*	2 5
*	2 8
*	3 5
*	3 8
*	5 8
*
*	Pass 2:
*
*	PVeMatches:
*	1 3 8 0 0 0 0 0
**	0 2 3 5 8 0 0 0
*	1 3 5 8 0 0 0 0
*	0 1 2 5 8 0 0 0
**	1 2 3 8 0 0 0 0
*	0 1 2 3 5 0 0 0
*
*	Updated PVe (Sorted):
*	0 1
*	0 2
*	0 3
*	0 5
*	0 8
**	1 2
*	1 3
*	1 5
*	1 8
*	2 3
*	2 5
*	2 8
*	3 5
*	3 8
*	5 8
**
*	Pass 3:
*
*	PVeMatches:
*	1 2 3 5 8 0 0 0
*	0 2 3 5 8 0 0 0
*	0 1 3 5 8 0 0 0
*	0 1 2 5 8 0 0 0
*	0 1 2 3 8 0 0 0
*	0 1 2 3 5 0 0 0
*
*	Updated PVe (Sorted):
*	No change => Done!
*
*
*	Code:
*
*	int *PVeTest, pvetest;
*
*	PVeTest = malloc(6*6*2 * sizeof *PVeTest);
**
*	PVeTest[0] = 1;  PVeTest[1] = 0;
*	PVeTest[2] = 2;  PVeTest[3] = 8;
*	PVeTest[4] = 1;  PVeTest[5] = 3;
*	PVeTest[6] = 8;  PVeTest[7] = 3;
*	PVeTest[8] = 3;  PVeTest[9] = 2;
*	PVeTest[10] = 8; PVeTest[11] = 1;
*	PVeTest[12] = 5; PVeTest[13] = 8;
*
*	pvetest = 7;
*
*	pve = pvetest;
*	array_print_i(pve,2,PVeTest);
**
*	FindPeriodicConnections(PVeTest,&pvetest,NVe); 
*
*	pve = pvetest;
*	array_print_i(pve,2,PVeTest);
*
*/

	void find_periodic_connections(int *PVe, int *pvePointer, int VeMax) {
	int i, j, k;
	int pve, *IndicesDummy, Modified, *PVe1D, NUnique, *PVeUnique, *PVeUniqueOver, *PVeMatches, *IndPVeMatches;
  int IndRow, LenF, match, row, col, n2, *PVePotential, IndPVe;

  pve = *pvePointer;

  //array_print_i(pve,2,PVe);

  // Sort Rows
  for (i = 0; i < pve; i++) PetscSortInt(2,&PVe[i*2+0]);

  // Sort Columns
  IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
  array_sort_i(pve,2,PVe,IndicesDummy,'R','T');
  free(IndicesDummy);

  //array_print_i(pve,2,PVe);

  Modified = 1;
  while (Modified) {
    Modified = 0;

    PVe1D = malloc(pve*2 * sizeof *PVe1D); // free
    k = 0;
    for (i = 0; i < pve; i++) {
    for (j = 0; j < 2; j++) {
      PVe1D[k] = PVe[i*2+j];
      k++;
    }}

    PetscSortInt(k,PVe1D);
    //array_print_i(1,k,PVe1D);
    
    NUnique = 1;
    PVeUniqueOver = malloc(pve*2 * sizeof *PVeUniqueOver); // free
    PVeUniqueOver[0] = PVe1D[0];
    for (i = 1; i < pve*2; i++) {
      if (PVe1D[i] != PVeUniqueOver[NUnique-1]) {
        PVeUniqueOver[NUnique] = PVe1D[i];
        NUnique++;
      }
    }
    free(PVe1D);

    PVeUnique = malloc(NUnique * sizeof *PVeUnique); // free
    for (i = 0; i < NUnique; i++) PVeUnique[i] = PVeUniqueOver[i];
    free(PVeUniqueOver);

    //array_print_i(1,NUnique,PVeUnique);

    PVeMatches    = malloc(NUnique*8 * sizeof *PVeMatches); // free
    IndPVeMatches = malloc(NUnique   * sizeof *IndPVeMatches); // free
    for (i = 0; i < NUnique*8; i++) PVeMatches[i]    = VeMax;
    for (i = 0; i < NUnique; i++)   IndPVeMatches[i] = 0;

    for (k = 0; k < 2; k++) {
      for (i = 0; i < pve; i++) {
        array_find_indexo_i(NUnique,PVeUnique,PVe[i*2+k],&IndRow,&LenF);

        for (j = 0, match = 0; j < IndPVeMatches[IndRow]; j++) {
          if (PVe[i*2+1-k] == PVeMatches[IndRow*8+j]) match = 1;
        }

        if (match == 0) {
          PVeMatches[IndRow*8+IndPVeMatches[IndRow]] = PVe[i*2+1-k];
       
          IndPVeMatches[IndRow]++;
          if (IndPVeMatches[IndRow] > 7) printf("Error: Too many periodic connections\n"), exit(1);
        }
      }
    }

    // sort rows
    for (i = 0; i < NUnique; i++) PetscSortInt(IndPVeMatches[i]+1,&PVeMatches[i*8+0]);

    //array_print_i(NUnique,8,PVeMatches);
    //array_print_i(1,NUnique,IndPVeMatches);

    PVePotential = malloc(2 * sizeof *PVePotential); // free

    for (row = 0; row < NUnique; row++) {
    for (col = 0; col < IndPVeMatches[row]; col++) {
      array_find_indexo_i(NUnique,PVeUnique,PVeMatches[row*8+col],&IndRow,&LenF);
      for (i = 0; i < IndPVeMatches[IndRow]; i++) {
        n2 = PVeMatches[IndRow*8+i];
        if (PVeUnique[row] != n2) {
          PVePotential[0] = PVeUnique[row];
          PVePotential[1] = n2;
          PetscSortInt(2,PVePotential);

          for (IndPVe = 0, match = 0; IndPVe < pve; IndPVe++) {
            if (PVe[IndPVe*2+0] == PVePotential[0] &&
                PVe[IndPVe*2+1] == PVePotential[1]) {
              match = 1;
              break;
            }
          }
         
          if (match == 0) {
            for (j = 0; j < 2; j++) PVe[pve*2+j] = PVePotential[j];
            pve++;
            Modified = 1;
          }
        }
      }
    }}
    free(PVeUnique);
    free(PVePotential);

    IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
    array_sort_i(pve,2,PVe,IndicesDummy,'R','T');
    free(IndicesDummy);

    //array_print_i(pve,2,PVe);

    free(PVeMatches);
    free(IndPVeMatches);
  }

  *pvePointer = pve;
}
