#include <stdlib.h>
#include <stdio.h>
//#include <string.h>

/*
 *	Purpose:
 *		Print arrays for debugging purposes.
 *
 *	Comments:
 *		If an elegant way is found to print arbitrary array types, revive array_print.
 *
 *	Notation:
 *
 *	References:
 *
 */

/*
void array_print(int m, int n, void *A, char *type)
{
	int i, j;

	long double *Ald = NULL;
	double      *Ad = NULL;
	float       *Af = NULL;
	long long   *All = NULL;
	long        *Al = NULL;
	int         *Ai = NULL;
	char        *Ac = NULL;

	if      (strstr(type,"ld")) Ald = A;
	else if (strstr(type,"d"))  Ad = A;
	else if (strstr(type,"f"))  Af = A;
	else if (strstr(type,"ll")) All = A;
	else if (strstr(type,"l"))  Al = A;
	else if (strstr(type,"i"))  Ai = A;
	else if (strstr(type,"c"))  Ac = A;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			if      (strstr(type,"ld")) printf("% .3Le ", Ald[i*n+j]);
			else if (strstr(type,"d"))  printf("% .3e ",  Ad[i*n+j]);
			else if (strstr(type,"f"))  printf("% .3f ",  Af[i*n+j]);
			else if (strstr(type,"ll")) printf("% 12lld ",All[i*n+j]);
			else if (strstr(type,"l"))  printf("% 12ld ", Al[i*n+j]);
			else if (strstr(type,"i"))  printf("% 12d ",  Ai[i*n+j]);
			else if (strstr(type,"c"))  printf("%c ",   Ac[i*n+j]);
		}
		printf("\n");
	}
	printf("\n");
}
*/

void array_print_i(int m, int n, int *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
			printf("% 12d ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

void array_print_l(int m, int n, long *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
		  printf("% 12ld ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

void array_print_ll(int m, int n, long long *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
		  printf("% 12lld ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

void array_print_f(int m, int n, float *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
		  printf("% .3f ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

void array_print_d(int m, int n, double *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
		  printf("% .3e ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

void array_print_ld(int m, int n, long double *A)
{
	int i, j;

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)
		  printf("% .3Le ",A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}
