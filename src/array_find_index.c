#include <stdlib.h>
#include <stdio.h>

/*
 *	Purpose:
 *	Find the pointer to the first index whose array entries correspond to the value to be found.
 *
 *	Comments:
 *	The 'O' in the function name indicates that the input array is (O)rdered. In this case, a
 *	divide and conquer algorithm is used.
 *	Note the usage of integer division (i.e. a = b/c is analogous to a = floor(b/c))
 *
 *	Notation:
 *	A   : Input array to looked in.
 *	val : Value of A to be found.
 *
 *	Example (array_find_indexo_i):
 *
 *	A = [0 1 1 4 5 7 7 7 7 7]
 *	val = 7
 *
 *	Pass 1: IndL = 0; InR = 9; IndM = 4; IdxF[0] = 5;
 *	Pass 2: IndL = 5; InR = 9; IndM = 7; IdxF[0] = 7;  Found val.
 *
 *	Code:
 *
 *	int *Atest,Idxtest,Lentest;
 *
 *	Atest = malloc(10 *sizeof *Atest);
 *	Atest[0] = 0; Atest[1] = 1; Atest[2] = 1; Atest[3] = 4; Atest[4] = 5;
 *	Atest[5] = 7; Atest[6] = 7; Atest[7] = 7; Atest[8] = 7; Atest[9] = 7;
 *
 *	array_find_indexo_i(10,Atest,7,Idxtest,Lentest);
 *	printf("%d %d\n",Idxtest,Lentest);
 *
 */

void array_find_indexo_i(int LenA, int *A, int val, int *IdxF, int *LenF) {
  int IndL, IndR, IndM, Aval;

  *LenF = 0;

  IndL = 0;
  IndR = LenA-1;
  IndM = (IndL + IndR)/2;
  Aval = A[IndM];

  if (Aval == val) *LenF = 1;
  while(Aval != val) {
    if      (val < Aval) IndR = IndM-1;
    else if (val > Aval) IndL = IndM+1;
    IndM = (IndL + IndR)/2;
    Aval = A[IndM];
  
    //printf("Loop: %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*LenF);

    if (Aval == val) {
      *LenF = 1;
      break;
    }

    if (IndL == IndR) {
      printf("Loop break: %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*LenF);
      break;
    }
  }
  if (*LenF == 0) printf("Did not find any matches in array_find_indexo_i"), exit(1);

  // Find how many times the value is repeated and the first entry
  IndL = IndM; IndR = IndM;
  while(IndL > 0      && A[IndL-1] == Aval) IndL--;
  while(IndR < LenA-1 && A[IndR+1] == Aval) IndR++;

  *LenF = IndR-IndL+1;
  *IdxF = IndL;

  //ArrayPrinti(1,LenA,A);
  //printf("%d %d %d %d %d %d %d\n",IndL,IndR,IndM,val,Aval,*IdxF,*LenF);
  
}
