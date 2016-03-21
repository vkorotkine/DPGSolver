#include <stdlib.h>

#include "mkl.h"
#include "petscsys.h"

/*
 *	Purpose:
 *		Sort multidimensional arrays maintaining fixed columns, with decreasing priority from first to last row.
 *
 *	Comments:
 *		The first row is sorted with a slightly different routine as it does not have an associated Master row.
 *		mkl_simatcopy is used to transpose arrays (i.e. convert from row major to column major storage and vice versa).
 *
 *	Notation:
 *		A       : Input array to be sorted of size NRows x NCols.
 *		Indices : Indices of the rearrangement.
 *
 *		M       : (M)aster
 *		S       : (S)lave
 *
 *	Example (array_sort_i):
 *		Notation: The array, A, is listed above and the Indices array is listed below the dashed line.
 *
 *		Unsorted:
 *		0 1 1 2 3 0 0 1 2 1
 *		1 3 3 5 1 0 1 0 4 3
 *		3 1 4 3 0 1 2 1 2 5
 *		-------------------
 *		0 1 2 3 4 5 6 7 8 9
 *
 *
 *		First row:
 *		0 0 0 1 1 1 1 2 2 3
 *		1 1 0 3 3 0 3 5 4 1
 *		3 2 1 5 4 1 1 3 2 0
 *		-------------------
 *		0 6 5 9 2 7 1 3 8 4
 *
 *		Second row:
 *		0-block               1-block               2-block               3-block (no change)
 *		0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
 *		0 1 1 3 3 0 3 5 4 1   0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
 *		1 2 3 5 4 1 1 3 2 0   1 2 3 1 4 5 1 3 2 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0
 *		-------------------   -------------------   -------------------   -------------------
 *		5 6 0 9 2 7 1 3 8 4   5 6 0 7 2 9 1 3 8 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4
 *
 *		Third row:
 *		00-block (no ch.)     01-block (no ch.)     10-block (no ch.)     13-block
 *		0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
 *		0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 5 4 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
 *		1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 4 5 1 2 3 0   1 2 3 1 1 4 5 2 3 0
 *		-------------------   -------------------   -------------------   -------------------
 *		5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 2 9 1 8 3 4   5 6 0 7 1 2 9 8 3 4
 *
 *		24-block (no ch.)     25-block (no ch.)     31-block (no ch.)     Sorted!
 *		0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3   0 0 0 1 1 1 1 2 2 3
 *		0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1   0 1 1 0 3 3 3 4 5 1
 *		1 2 3 1 1 4 5 2 3 0   1 2 3 1 1 4 5 2 3 0   1 2 3 1 1 4 5 2 3 0
 *		-------------------   -------------------   -------------------
 *		5 6 0 7 1 2 9 8 3 4   5 6 0 7 1 2 9 8 3 4   5 6 0 7 1 2 9 8 3 4
 *
 *		Code:
 *
 *			int *Atest, *Indicestest, Rtest = 3, Ctest = 10;
 *
 *			Atest = malloc(Rtest*Ctest *sizeof *Atest);
 *			Atest[0]  = 0; Atest[1]  = 1; Atest[2]  = 1; Atest[3]  = 2; Atest[4]  = 3;
 *			Atest[5]  = 0; Atest[6]  = 0; Atest[7]  = 1; Atest[8]  = 2; Atest[9]  = 1;
 *			Atest[10] = 1; Atest[11] = 3; Atest[12] = 3; Atest[13] = 5; Atest[14] = 1;
 *			Atest[15] = 0; Atest[16] = 1; Atest[17] = 0; Atest[18] = 4; Atest[19] = 3;
 *			Atest[20] = 3; Atest[21] = 1; Atest[22] = 4; Atest[23] = 3; Atest[24] = 0;
 *			Atest[25] = 1; Atest[26] = 2; Atest[27] = 1; Atest[28] = 2; Atest[29] = 5;
 *
 *			Indicestest = malloc(Ctest *sizeof *Indicestest);
 *			for (i = 0; i < Ctest; i++) Indicestest[i] = i;
 *
 *			array_print_i(Rtest,Ctest,Atest);
 *			array_print_i(1,Ctest,Indicestest);
 *
 *			array_sort_i(Rtest,Ctest,Atest,Indicestest,'R','N');
 *
 *			array_print_i(Rtest,Ctest,Atest);
 *			array_print_i(1,Ctest,Indicestest);
 *
 *
 *		Code to test array_sort_d:
 *
 *			int *Indicestest, Rtest = 3, Ctest = 10;
 *			double *Atest;
 *
 *			Atest = malloc(Rtest*Ctest * sizeof *Atest);
 *			Atest[0] = 0.1; Atest[1] = 1.1; Atest[2] = 1.1; Atest[3] = 2.1; Atest[4] = 3.1;
 *			Atest[5] = 0.1; Atest[6] = 0.1; Atest[7] = 1.1; Atest[8] = 2.1; Atest[9] = 1.1;
 *			Atest[10] = 1.1; Atest[11] = 3.1; Atest[12] = 3.1; Atest[13] = 5.1; Atest[14] = 1.1;
 *			Atest[15] = 0.1; Atest[16] = 1.1; Atest[17] = 0.1; Atest[18] = 4.1; Atest[19] = 3.1;
 *			Atest[20] = 3.1; Atest[21] = 1.1; Atest[22] = 4.1; Atest[23] = 3.1; Atest[24] = 0.1;
 *			Atest[25] = 1.1; Atest[26] = 2.1; Atest[27] = 1.1; Atest[28] = 2.1; Atest[29] = 5.1;
 *
 *			Indicestest = malloc(Ctest * sizeof *Indicestest);
 *			for (i = 0; i < Ctest; i++) Indicestest[i] = i;
 *
 *			array_print_d(Rtest,Ctest,Atest);
 *			array_print_i(1,Ctest,Indicestest);
 *
 *			array_sort_d(Rtest,Ctest,Atest,Indicestest,'R','N');
 *
 *			array_print_d(Rtest,Ctest,Atest);
 *			array_print_i(1,Ctest,Indicestest);
 *
 *	References:
 *
 */

void array_sort_i(int NRows, int NCols, int *A, int *Indices, char ordering, char trans)
{
	int i, count, swap,
	    row, col, colMs[NRows], colMe[NRows], SortLen, rowM, rowS,
	    *Indicestmp, *IndicesInter, *Atmp;

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_simatcopy(ordering,trans,NRows,NCols,1.,(float *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	} else {
		printf("Error: Invalid ordering/trans input to array_sort_i"), exit(1);
	}
// array_print_i(NRows,NCols,A);

	for (row = 0; row < NRows; row++) {
		if (row == 0) {
			colMs[row] = 0;
			colMe[row] = NCols-1;
			SortLen = colMe[row] - colMs[row] + 1;

			IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
			Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
			for (i = 0; i < SortLen; i++)
				Indicestmp[i] = i;

			PetscSortIntWithArray(SortLen,&A[row*NCols],&Indicestmp[0]);
			for (i = 0; i < SortLen; i++) IndicesInter[i] = Indices[Indicestmp[i]];
			for (i = 0; i < SortLen; i++) Indices[i]      = IndicesInter[i];
			free(IndicesInter);

			Atmp = malloc(NCols * sizeof *Atmp); // free
			// Rearrange other rows accordingly
			for (rowM = row+1; rowM < NRows; rowM++) {
				for (col = 0; col < NCols; col++) Atmp[col]         = A[rowM*NCols+col];
				for (col = 0; col < NCols; col++) A[rowM*NCols+col] = Atmp[Indicestmp[col]];
			}
		free(Indicestmp);
		free(Atmp);
		} else {
			colMs[0] = colMe[row-1] = 0;
			while (colMe[row-1] < NCols-1) {
				for (rowM = 0; rowM < row; rowM++) {
					if (rowM > 0)
						colMs[rowM] = colMs[rowM-1];
					col = colMs[rowM]+1;
					while (col < NCols && A[rowM*NCols+col] == A[rowM*NCols+colMs[rowM]])
						col++;
					colMe[rowM] = col-1;
				}
				rowM -= 1;
				SortLen = colMe[rowM] - colMs[rowM]+1;

//printf("%d %d %d %d %d %d %d \n",row,rowM,colMs[0],colMe[0],colMs[rowM],colMe[rowM],SortLen);

				IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
				Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
				for (i = 0; i < SortLen; i++)
					Indicestmp[i] = i;

				PetscSortIntWithArray(SortLen,&A[row*NCols+colMs[rowM]],&Indicestmp[0]);
				for (i = 0; i < SortLen; i++) IndicesInter[i]        = Indices[Indicestmp[i]+colMs[rowM]];
				for (i = 0; i < SortLen; i++) Indices[i+colMs[rowM]] = IndicesInter[i];
				free(IndicesInter);

				Atmp = malloc(SortLen * sizeof *Atmp); // free

				// Rearrange slave rows accordingly
				for (rowS = row+1; rowS < NRows; rowS++) {
					for (count = 0, col = colMs[rowM]; col <= colMe[rowM]; col++) {
						Atmp[count]         = A[rowS*NCols+col];
						count++;
					}
					for (count = 0, col = colMs[rowM]; col <= colMe[rowM]; col++) {
						A[rowS*NCols+col] = Atmp[Indicestmp[count]];
						count++;
					}
				}
				free(Indicestmp);
				free(Atmp);

				colMs[0] = colMe[row-1]+1;
			}
		}
// array_print_i(NRows,NCols,A);
// array_print_i(1,NCols,Indices);
	}

// array_print_i(NRows,NCols,A);

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_simatcopy(ordering,trans,NRows,NCols,1.,(float *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	}
// array_print_i(NRows,NCols,A);
}



void array_sort_d(int NRows, int NCols, double *A, int *Indices, char ordering, char trans)
{
	int    i, count, swap,
	       row, col, colMs[NRows], colMe[NRows], SortLen, rowM, rowS,
	       *Indicestmp, *IndicesInter;
	double *Atmp;

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_dimatcopy(ordering,trans,NRows,NCols,1.,(double *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	} else {
		printf("Error: Invalid ordering/trans input to array_sort_d"), exit(1);
	}
// array_print_d(NRows,NCols,A);

	for (row = 0; row < NRows; row++) {
		if (row == 0) {
			colMs[row] = 0;
			colMe[row] = NCols-1;
			SortLen = colMe[row] - colMs[row] + 1;

			IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
			Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
			for (i = 0; i < SortLen; i++)
				Indicestmp[i] = i;

			PetscSortRealWithPermutation(SortLen,&A[row*NCols],&Indicestmp[0]);
			for (i = 0; i < SortLen; i++) IndicesInter[i]  = Indices[Indicestmp[i]];
			for (i = 0; i < SortLen; i++) Indices[i]       = IndicesInter[i];
			free(IndicesInter);

			Atmp = malloc(NCols * sizeof *Atmp); // free
			// Rearrange rows (including the current row) accordingly
			for (rowM = row; rowM < NRows; rowM++) {
				for (col = 0; col < NCols; col++) Atmp[col]         = A[rowM*NCols+col];
				for (col = 0; col < NCols; col++) A[rowM*NCols+col] = Atmp[Indicestmp[col]];
			}
			free(Indicestmp);
			free(Atmp);
		} else {
			colMs[0] = colMe[row-1] = 0;
			while (colMe[row-1] < NCols-1) {
				for (rowM = 0; rowM < row; rowM++) {
					if (rowM > 0)
						colMs[rowM] = colMs[rowM-1];
					col = colMs[rowM]+1;
					while (col < NCols && A[rowM*NCols+col] == A[rowM*NCols+colMs[rowM]])
						col++;
					colMe[rowM] = col-1;
				}
				rowM -= 1;
				SortLen = colMe[rowM] - colMs[rowM]+1;

// printf("%d %d %d %d %d %d %d \n",row,rowM,colMs[0],colMe[0],colMs[rowM],colMe[rowM],SortLen);

				IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
				Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
				for (i = 0; i < SortLen; i++)
					Indicestmp[i] = i;

				PetscSortRealWithPermutation(SortLen,&A[row*NCols+colMs[rowM]],&Indicestmp[0]);
				for (i = 0; i < SortLen; i++) IndicesInter[i]        = Indices[Indicestmp[i]+colMs[rowM]];
				for (i = 0; i < SortLen; i++) Indices[i+colMs[rowM]] = IndicesInter[i];
				free(IndicesInter);

				Atmp = malloc(SortLen * sizeof *Atmp); // free

				// Rearrange rows (including the current row) accordingly
				for (rowS = row; rowS < NRows; rowS++) {
					for (count = 0, col = colMs[rowM]; col <= colMe[rowM]; col++) {
						Atmp[count]         = A[rowS*NCols+col];
						count++;
					}
					for (count = 0, col = colMs[rowM]; col <= colMe[rowM]; col++) {
						A[rowS*NCols+col] = Atmp[Indicestmp[count]];
						count++;
					}
				}
				free(Indicestmp);
				free(Atmp);

				colMs[0] = colMe[row-1]+1;
			}
		}
// printf("row: %d: \n",row);
// array_print_d(NRows,NCols,A);
// array_print_i(1,NCols,Indices);
	}

// array_print_d(NRows,NCols,A);

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_dimatcopy(ordering,trans,NRows,NCols,1.,(double *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	}
// array_print_d(NRows,NCols,A);
}
