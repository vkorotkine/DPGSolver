// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "parameters.h"

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
			if (fabs(A[i]) > norm) norm = fabs(A[i]);
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]);
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
		if (isnan(A[i]))
			printf("Error: Entry in array is 'nan'.\n"), EXIT_MSG;
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
			if (fabs(As[i]-Bs[i]) > norm_num) norm_num = fabs(As[i]-Bs[i]);
			if (fabs(As[i])       > norm_den) norm_den = fabs(As[i]);
		}
	} else if (strstr(NormType,"L1")) {
		for (i = 0; i < LenA; i++) {
			norm_num += fabs(As[i]-Bs[i]);
			norm_den += fabs(As[i]);
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
