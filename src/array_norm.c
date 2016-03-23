#include <string.h>

/*
 *	Purpose:
 *	Compute various norms.
 *
 *	Comments:
 *
 *	Notation:
 *
 *  References:
 *
 */

double array_norm_d(int LenA, double *A, char *NormType)
{
	int    i;
	double norm;

	if (strstr(NormType,"Inf") != NULL) {
		norm = 0.;
		for (i = 0; i < LenA; i++)
			if (fabs(A[i]) > norm) norm = fabs(A[i]);
	} else if (strstr(NormType,"L1") != NULL) {
	} else if (strstr(NormType,"L2") != NULL) {
	}

	return norm;
}
