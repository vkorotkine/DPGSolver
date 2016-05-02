#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return the "matrix" ChiRef_rst representing the orthonormal basis functions evaluated at the provided quadrature
 *		nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrix" is returned as a 1D array.
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified (Possibly Bergot, Chan)
 */

double *basis_PYR(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                  const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, Indbf, Indn, Nbf, mu_ij, IndChi;
	double       *ChiRef_rst, *a, *b, *c, a_n, b_n, c_n;

	if (d != 3)
		printf("Error: basis_PYR only supports d = 3.\n"), exit(1);

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_PYR(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (1.0/6.0*((P+1)*(P+2)*(2*P+3)));

	ChiRef_rst = malloc(Nn*Nbf * sizeof *ChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P; i <= iMax; i++) {
	for (j = 0, jMax = P; j <= jMax; j++) {
		mu_ij = max(i,j);
		for (k = 0, kMax = P-mu_ij; k <= kMax; k++) {
			printf("%d %d %d %d %d %d\n",Indbf,i,j,k,i+j+k,mu_ij+k);
			for (Indn = 0; Indn < Nn; Indn++) {
				IndChi = Indbf*Nn+Indn;
				a_n = a[Indn];
				b_n = b[Indn];
				c_n = c[Indn];

// Orthonormal (WV nodes of higher order than expected)
				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) mu_ij);
				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

// Orthonormal (GL nodes)
//				ChiRef_rst[IndChi]  = pow(1.0-c_n,(double) mu_ij);
//				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

// Orthonormal (WV nodes)
//				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) (i+j));
//				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

// Orthonormal (GL nodes)
//				ChiRef_rst[IndChi]  = pow(1.0-c_n,(double) (i+j));
//				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

// Remaining part of the basis
				ChiRef_rst[IndChi]  = 1.0;
				ChiRef_rst[IndChi] *= jacobiP(a_n,0.0,0.0,(double) i);
				ChiRef_rst[IndChi] *= jacobiP(b_n,0.0,0.0,(double) j);

/*
if (Indbf < 10) {
				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) mu_ij);
				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);
//				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) (i+j));
//				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);
				ChiRef_rst[IndChi] *= jacobiP(a_n,0.0,0.0,(double) i);
				ChiRef_rst[IndChi] *= jacobiP(b_n,0.0,0.0,(double) j);
} else if (Indbf == 10) {
				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) 2);
				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(2+1.0),0.0,(double) 0);
				ChiRef_rst[IndChi] *= jacobiP(a_n,0.0,0.0,(double) 2);
				ChiRef_rst[IndChi] *= jacobiP(b_n,0.0,0.0,(double) 2);
} else if (Indbf == 11) {
	ChiRef_rst[IndChi] = 1.0/4.0*(3.0*pow(a_n,2.0)-1)*(3.0*pow(b_n,2.0)-1.0)*pow(1.0-c_n,2.0);
} else if (Indbf == 12) {
	ChiRef_rst[IndChi] = 1.0/4.0*(3.0*pow(a_n,2.0)-1)*(3.0*pow(b_n,2.0)-1.0)*pow(1.0-c_n,4.0);
} else if (Indbf == 13) {
				ChiRef_rst[IndChi]  = pow(2.0,1.25)*pow(1.0-c_n,(double) (2+2));
				ChiRef_rst[IndChi] *= jacobiP(c_n,2.0*(2+2+1.0),0.0,(double) 0);
				ChiRef_rst[IndChi] *= jacobiP(a_n,0.0,0.0,(double) 2);
				ChiRef_rst[IndChi] *= jacobiP(b_n,0.0,0.0,(double) 2);
} else {
	ChiRef_rst[IndChi] = 0.0;
}
*/

		}
		Indbf++;
	}}}

	// Transpose ChiRef_rst
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiRef_rst,Nn,Nbf);

//	array_print_d(Nn,Nbf,ChiRef_rst,'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return ChiRef_rst;
}
