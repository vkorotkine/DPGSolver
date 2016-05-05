#include <stdlib.h>

#include "functions.h"

#include "mkl.h"
#include "petscsys.h"

/*
 *	Purpose:
 *		Sort multidimensional arrays maintaining fixed columns, with decreasing priority from first to last row.
 *
 *	Comments:
 *		The first row is sorted with a slightly different routine as it does not have an associated Master row.
 *		mkl_simatcopy is used to transpose arrays (i.e. convert from row major to column major storage and vice versa).
 *		When finished coding, check if array_sort_i is used; if not, delete it. (ToBeDeleted)
 *
 *	Notation:
 *		A       : Input array to be sorted of size NRows x NCols.
 *		Indices : Indices of the rearrangement.
 *
 *		M       : (M)aster
 *		S       : (S)lave
 *
 *	References:
 *
 */

void array_sort_ui(unsigned int NRows, unsigned int NCols, unsigned int *A, unsigned int *Indices, const char ordering,
                   const char trans)
{
	unsigned int i, count, swap,
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
// array_print_i(NRows,NCols,A,'R');

	for (row = 0; row < NRows; row++) {
		if (row == 0) {
			colMs[row] = 0;
			colMe[row] = NCols-1;
			SortLen = colMe[row] - colMs[row] + 1;

			IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
			Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
			for (i = 0; i < SortLen; i++)
				Indicestmp[i] = i;

			PetscSortIntWithArray(SortLen,(int *)&A[row*NCols],(int *)Indicestmp);
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
					while (col < NCols && A[rowM*NCols+col] == A[rowM*NCols+colMs[rowM]]) {
						// Ensure that all rows are taken into account and not only the previous row
						if (rowM > 0 && col > colMe[rowM-1])
							break;
						col++;
					}
					colMe[rowM] = col-1;
				}
				rowM -= 1;
				SortLen = colMe[rowM] - colMs[rowM]+1;

//printf("%d %d %d %d %d %d %d %d \n",row,col,rowM,colMs[0],colMe[0],colMs[rowM],colMe[rowM],SortLen);

				IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
				Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
				for (i = 0; i < SortLen; i++)
					Indicestmp[i] = i;

				PetscSortIntWithArray(SortLen,(int *)&A[row*NCols+colMs[rowM]],(int *)Indicestmp);
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
//array_print_ui(NCols,NRows,A,'C');
//array_print_ui(1,NCols,Indices,'R');
	}

// array_print_i(NRows,NCols,A,'R');

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_simatcopy(ordering,trans,NRows,NCols,1.,(float *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	}
// array_print_i(NRows,NCols,A,'R');
}



void array_sort_i(unsigned int NRows, unsigned int NCols, int *A, unsigned int *Indices, const char ordering,
                  const char trans)
{
	unsigned int i, count, swap,
	             row, col, colMs[NRows], colMe[NRows], SortLen, rowM, rowS,
	             *Indicestmp, *IndicesInter;
	int          *Atmp;

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_simatcopy(ordering,trans,NRows,NCols,1.,(float *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	} else {
		printf("Error: Invalid ordering/trans input to array_sort_i"), exit(1);
	}
// array_print_i(NRows,NCols,A,'R');

	for (row = 0; row < NRows; row++) {
		if (row == 0) {
			colMs[row] = 0;
			colMe[row] = NCols-1;
			SortLen = colMe[row] - colMs[row] + 1;

			IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
			Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
			for (i = 0; i < SortLen; i++)
				Indicestmp[i] = i;

			PetscSortIntWithArray(SortLen,&A[row*NCols],(int *)Indicestmp);
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
					while (col < NCols && A[rowM*NCols+col] == A[rowM*NCols+colMs[rowM]]) {
						// Ensure that all rows are taken into account and not only the previous row
						if (rowM > 0 && col > colMe[rowM-1])
							break;
						col++;
					}
					colMe[rowM] = col-1;
				}
				rowM -= 1;
				SortLen = colMe[rowM] - colMs[rowM]+1;

//printf("%d %d %d %d %d %d %d \n",row,rowM,colMs[0],colMe[0],colMs[rowM],colMe[rowM],SortLen);

				IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
				Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
				for (i = 0; i < SortLen; i++)
					Indicestmp[i] = i;

				PetscSortIntWithArray(SortLen,&A[row*NCols+colMs[rowM]],(int *)Indicestmp);
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
// array_print_i(NRows,NCols,A,'R');
// array_print_i(1,NCols,Indices,'R');
	}

// array_print_i(NRows,NCols,A,'R');

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_simatcopy(ordering,trans,NRows,NCols,1.,(float *) A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	}
// array_print_i(NRows,NCols,A,'R');
}



void array_sort_d(unsigned int NRows, unsigned int NCols, double *A, unsigned int *Indices, const char ordering,
                  const char trans)
{
	unsigned int i, count, swap,
	             row, col, colMs[NRows], colMe[NRows], SortLen, rowM, rowS,
	             *Indicestmp, *IndicesInter;
	double       *Atmp;

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_dimatcopy(ordering,trans,NRows,NCols,1.,A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	} else {
		printf("Error: Invalid ordering/trans input to array_sort_d"), exit(1);
	}
// array_print_d(NRows,NCols,A,'R');

	for (row = 0; row < NRows; row++) {
		if (row == 0) {
			colMs[row] = 0;
			colMe[row] = NCols-1;
			SortLen = colMe[row] - colMs[row] + 1;

			IndicesInter = malloc(SortLen * sizeof *IndicesInter); // free
			Indicestmp   = malloc(SortLen * sizeof *Indicestmp); // free
			for (i = 0; i < SortLen; i++)
				Indicestmp[i] = i;

			PetscSortRealWithPermutation(SortLen,&A[row*NCols],(int *)Indicestmp);
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

				PetscSortRealWithPermutation(SortLen,&A[row*NCols+colMs[rowM]],(int *)Indicestmp);
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
// array_print_d(NRows,NCols,A,'R');
// array_print_i(1,NCols,Indices,'R');
	}

// array_print_d(NRows,NCols,A,'R');

	if ((ordering == 'R' && trans == 'T') || (ordering == 'C' && trans == 'N')) {
		mkl_dimatcopy(ordering,trans,NRows,NCols,1.,A,NCols,NRows);
		swap = NRows; NRows = NCols; NCols = swap;
	} else if ((ordering == 'R' && trans == 'N') || (ordering == 'C' && trans == 'T')) {
		// Don't do anything.
	}
// array_print_d(NRows,NCols,A,'R');
}
