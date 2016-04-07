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
		norm = 0.;
		for (i = 0; i < LenA; i++)
			norm += fabs(A[i]);
	} else if (strstr(NormType,"L2") != NULL) {
		norm = 0.;
		for (i = 0; i < LenA; i++)
			norm += pow(A[i],2);
		norm = sqrt(norm);
	}

	return norm;
}
