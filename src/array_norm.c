#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *	Purpose:
 *		Compute various norms.
 *
 *	Comments:
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

	if (strstr(NormType,"Inf") != NULL) {
		for (i = 0; i < LenA; i++)
			if (fabs(A[i]) > norm) norm = fabs(A[i]);
	} else if (strstr(NormType,"L1") != NULL) {
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]);
	} else if (strstr(NormType,"L2") != NULL) {
		printf("Error: L2 norm not supported for unsigned int (norm).\n"), exit(1);
	}

	return norm;
}

double array_norm_d(const unsigned int LenA, const double *A, const char *NormType)
{
	unsigned int i;
	double       norm = 0.0;

	if (strstr(NormType,"Inf") != NULL) {
		for (i = 0; i < LenA; i++)
			if (fabs(A[i]) > norm) norm = fabs(A[i]);
	} else if (strstr(NormType,"L1") != NULL) {
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]);
	} else if (strstr(NormType,"L2") != NULL) {
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
	int *As = (int *) A,
	    *Bs = (int *) B;
	int norm = 0;

	if (strstr(NormType,"Inf") != NULL) {
		for (i = 0; i < LenA; i++)
			if (fabs(As[i]-Bs[i]) > norm) norm = fabs(As[i]-Bs[i]);
	} else if (strstr(NormType,"L1") != NULL) {
		for (i = 0; i < LenA; i++)
			norm += fabs(As[i]-Bs[i]);
	} else if (strstr(NormType,"L2") != NULL) {
		printf("Error: L2 norm not supported for unsigned int (norm diff).\n"), exit(1);
	}

	return (unsigned int) norm;
}

double array_norm_diff_d(const unsigned int LenA, const double *A, const double *B, const char *NormType)
{
	unsigned int i;
	double       norm = 0.0;

	if (strstr(NormType,"Inf") != NULL) {
		for (i = 0; i < LenA; i++)
			if (fabs(A[i]-B[i]) > norm) norm = fabs(A[i]-B[i]);
	} else if (strstr(NormType,"L1") != NULL) {
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]-B[i]);
	} else if (strstr(NormType,"L2") != NULL) {
		for (i = 0; i < LenA; i++)
			norm += pow(A[i]-B[i],2);
		norm = sqrt(norm);
	}

	return norm;
}
