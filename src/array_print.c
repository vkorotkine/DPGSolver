#include <stdlib.h>
#include <stdio.h>
//#include <string.h>

/*
 *	Purpose:
 *		Print arrays for debugging purposes.
 *
 *	Comments:
 *		If an elegant way is found to print arbitrary array types, revive array_print (ToBeDeleted).
 *
 *	Notation:
 *
 *	References:
 *
 */

/*
void array_print(const unsigned int m, const unsigned int n, void *A, char *type)
{
	unsigned int i, j;

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

void array_print_ui(const unsigned int m, const unsigned int n, const unsigned int *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_i(const unsigned int m, const unsigned int n, int *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12d ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_l(const unsigned int m, const unsigned int n, long *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12ld ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12ld ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_ll(const unsigned int m, const unsigned int n, long long *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12lld ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% 12lld ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_f(const unsigned int m, const unsigned int n, float *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3f ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3f ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_d(const unsigned int m, const unsigned int n, double *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3e ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3e ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}

void array_print_ld(const unsigned int m, const unsigned int n, long double *A, const char layout)
{
	unsigned int i, j;

	switch (layout) {
	case 'R':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3Le ",A[i*n+j]);
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++)
				printf("% .3Le ",A[i+j*m]);
			printf("\n");
		}
		printf("\n");
		break;
	}
}
