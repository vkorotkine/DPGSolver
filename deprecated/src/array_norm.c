// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "array_norm.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "petscmat.h"

#include "Macros.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute various norms.
 *
 *	Comments:
 *		array_norm_diff_d computes the relative error if round-off does not affect the result.
 *
 *	Notation:
 *
 *	References:
 *
 */

unsigned int array_norm_ui(const unsigned int LenA, const unsigned int *A, const char *NormType)
{
	unsigned int i;
	unsigned int norm = 0;

	if (strstr(NormType,"Inf")) {
		for (i = 0; i < LenA; i++)
			if (A[i] > norm) norm = A[i];
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++)
			norm += A[i];
	} else if (strstr(NormType,"L2")) {
		printf("Error: L2 norm not supported for unsigned int (norm).\n"), exit(1);
	}

	return norm;
}

double array_norm_d(const unsigned int LenA, const double *A, const char *NormType)
{
	unsigned int i;
	double       norm = 0.0;

	for (i = 0; i < LenA; i++) {
		if (isnan(A[i])) {
			array_print_d(1,LenA,A,'R');
			printf("Error: Entry in array is 'nan'.\n"), EXIT_MSG;
		}
	}

	if (strstr(NormType,"Inf")) {
		for (i = 0; i < LenA; i++)
			if (fabs(A[i]) > norm) norm = fabs(A[i]);
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]);
	} else if (strstr(NormType,"L2")) {
		for (i = 0; i < LenA; i++)
			norm += pow(A[i],2);
		norm = sqrt(norm);
	}

	return norm;
}

unsigned int array_norm_diff_ui(const unsigned int LenA, const unsigned int *A, const unsigned int *B,
                                const char *NormType)
{
	unsigned int i;
	int          *As = (int *) A,
	             *Bs = (int *) B,
	             norm_num = 0, norm_den = 0;

	if (strstr(NormType,"Inf")) {
		for (i = 0; i < LenA; i++) {
			if (abs(As[i]-Bs[i]) > norm_num) norm_num = abs(As[i]-Bs[i]);
			if (abs(As[i])       > norm_den) norm_den = abs(As[i]);
		}
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++) {
			norm_num += abs(As[i]-Bs[i]);
			norm_den += abs(As[i]);
		}
	} else if (strstr(NormType,"L2")) {
		printf("Error: L2 norm not supported.\n"), EXIT_MSG;
	}

	return (unsigned int) norm_num;
}

double array_norm_diff_d(const unsigned int LenA, const double *A, const double *B, const char *NormType)
{
	unsigned int i;
	double       norm_num = 0.0, norm_den = 0.0;

	for (i = 0; i < LenA; i++) {
		if (isnan(A[i]) || isnan(B[i]))
			printf("Error: Entry in array is 'nan'.\n"), EXIT_MSG;
	}

	if (strstr(NormType,"Inf")) {
		for (i = 0; i < LenA; i++) {
			if (fabs(A[i]-B[i]) > norm_num) norm_num = fabs(A[i]-B[i]);
			if (fabs(A[i])      > norm_den) norm_den = fabs(A[i]);
		}
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++) {
			norm_num += fabs(A[i]-B[i]);
			norm_den += fabs(A[i]);
		}
	} else if (strstr(NormType,"L2")) {
		for (i = 0; i < LenA; i++) {
			norm_num += pow(A[i]-B[i],2);
			norm_den += pow(A[i],2);
		}
		norm_num = sqrt(norm_num);
		norm_den = sqrt(norm_den);
	}

	if (norm_den > 1.0) {
		return norm_num/norm_den;
	}
	else
		return norm_num;
}

double array_norm_diff_dc(const unsigned int LenA, const double *A, const double complex *B, const char *NormType)
{
	double *Br = malloc(LenA * sizeof *Br); // free

	for (size_t i = 0; i < LenA; i++)
		Br[i] = creal(B[i]);

	double norm_diff = array_norm_diff_d(LenA,A,Br,NormType);
	free(Br);

	return norm_diff;
}

static double array_norm_diff_d_sp (const unsigned int LenA, const double*const A, const double*const B,
                                    const PetscInt*const colsA, const PetscInt*const colsB)
{
	/*
	 *	Comments:
	 *		LenA must be less than or equal to LenB.
	 *		Currently only supported the relative infinity norm.
	 */

	double norm_num = 0.0,
	       norm_den = 0.0;
	for (size_t i = 0, j = 0; i < LenA; i++) {
		while (colsA[i] != colsB[j])
			j++;

		double tmp = fabs(A[i]-B[j]);
		if (tmp > norm_num)
			norm_num = tmp;

		tmp = fabs(A[i]);
		if (tmp > norm_den)
			norm_den = tmp;
	}

	return norm_num/norm_den;
}

double PetscMatAIJ_norm_diff_d(const unsigned int NRows, Mat A, Mat B, const char *NormType, const bool allow_uninit)
{
	/*
	 *	Comments:
	 *		It is possible that some entries of one of the two input matrices are not initialized when separating the
	 *		VOLUME and FACE LHS contributions to A. The (allow)_(un)(init)ialized flag, causes the norm to ignore
	 *		uninitialized entries if present.
	 */

	unsigned int i;
	double       norm_row, norm;

	int               ncols[2];
	const PetscInt    *cols[2];
	const PetscScalar *vals[2];

	norm = 0.0;
	if (strstr(NormType,"Inf")) {
		for (i = 0; i < NRows; i++) {

			MatGetRow(A,i,&ncols[0],&cols[0],&vals[0]);
			MatGetRow(B,i,&ncols[1],&cols[1],&vals[1]);

			if (!allow_uninit) {
				if (ncols[0] != ncols[1]) {
					printf("Error: Different number of non-zero columns in A (%d) and B (%d) on line %d.\n",
				           ncols[0],ncols[1],i);
					EXIT_UNSUPPORTED;
				}

				norm_row = array_norm_diff_d(ncols[0],vals[0],vals[1],"Inf");
			} else {
				const size_t indA = ( ncols[0] < ncols[1] ? 0 : 1 ),
				             indB = ( indA ? 0 : 1 );
				norm_row = array_norm_diff_d_sp(ncols[indA],vals[indA],vals[indB],cols[indA],cols[indB]);
			}

			if (norm_row > norm)
				norm = norm_row;

			MatRestoreRow(A,i,&ncols[0],&cols[0],&vals[0]);
			MatRestoreRow(A,i,&ncols[1],&cols[1],&vals[1]);
		}
	} else {
		printf("Error: Only infinity norm is supported.\n"), EXIT_MSG;
	}

	return norm;
}
