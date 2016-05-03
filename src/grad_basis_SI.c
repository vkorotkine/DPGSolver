#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return the "matrices" GradChiRef_rst representing the gradients of the orthonormal basis functions evaluated at
 *		the provided quadrature nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrices" are returned as a 1D array.
 *
 *	Notation:
 *
 *	References:
 *		ToBeModified: Add reference to thesis where it is explained how these derivatives are computed.
 */

double **grad_basis_SI(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                       const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, pder, Indbf, Indn, IndGrad, Nbf;
	double       **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, jPa, jPb, jPc, djPa, djPb, djPc;
	double       GradChiRef_tmp[3][5];

	if (d < 2 || d > 3)
		printf("Error: grad_basis_SI only supports d = [2,3].\n"), exit(1);

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_SI(Nn,d,rst,a,b,c);

	Nbf = factorial_ull(d+P)/(factorial_ull(d)*factorial_ull(P));

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		 GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P;   i <= iMax; i++) {
	for (j = 0, jMax = P-i; j <= jMax; j++) {
	for (k = 0, kMax = (d-2)*(P-i-j); k <= kMax; k++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			IndGrad = Indbf*Nn+Indn;
			a_n = a[Indn];
			b_n = b[Indn];
			c_n = c[Indn];

			jPa = jacobiP(a_n,0.0,          0.0,(double) i);
			jPb = jacobiP(b_n,2.0*i+1.0,    0.0,(double) j);
			jPc = jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

			djPa = grad_jacobiP(a_n,0.0,          0.0,(double) i);
			djPb = grad_jacobiP(b_n,2.0*i+1.0,    0.0,(double) j);
			djPc = grad_jacobiP(c_n,2.0*(i+j+1.0),0.0,(double) k);

			for (dim = 0; dim < d; dim++)
				for (pder = 0; pder < 5; pder++ )
					GradChiRef_tmp[dim][pder] = 0.0;

			// Obtain contributions from each partial derivative
			GradChiRef_tmp[0][0] = djPa* jPb;
			GradChiRef_tmp[1][0] = djPa* jPb;
			GradChiRef_tmp[1][1] =  jPa*djPb;
			GradChiRef_tmp[1][2] =  jPa* jPb;

			if (d == 3) {
				// djPa = 0 when i = 0
				// djPb = 0 when j = 0
				// djPc = 0 when k = 0
				GradChiRef_tmp[0][0] *= jPc;
				GradChiRef_tmp[1][0] *= jPc;
				GradChiRef_tmp[1][1] *= jPc;
				GradChiRef_tmp[1][2] *= jPc;
				GradChiRef_tmp[2][0]  = djPa* jPb* jPc;
				GradChiRef_tmp[2][1]  =  jPa*djPb* jPc;
				GradChiRef_tmp[2][2]  =  jPa* jPb* jPc;
				GradChiRef_tmp[2][3]  =  jPa* jPb*djPc;
				GradChiRef_tmp[2][4]  =  jPa* jPb* jPc;
			}

			GradChiRef_tmp[0][0] *=  2.0                  ;
			GradChiRef_tmp[1][0] *=  2.0/3.0*sqrt(3.0)*a_n;
			GradChiRef_tmp[1][1] *=  2.0/3.0*sqrt(3.0)    ;
			GradChiRef_tmp[1][2] *= -2.0/3.0*sqrt(3.0)*i  ;
			if (i > 0) {
				GradChiRef_tmp[0][0] *= pow(1.0-b_n,i-1.0);
				GradChiRef_tmp[1][0] *= pow(1.0-b_n,i-1.0);
				GradChiRef_tmp[1][2] *= pow(1.0-b_n,i-1.0);
			} else {
				GradChiRef_tmp[0][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
				GradChiRef_tmp[1][0] *= 0.0; // redundant (i = 0 -> djPa = 0)
				GradChiRef_tmp[1][2] *= 0.0;
			}
			GradChiRef_tmp[1][1] *= pow(1.0-b_n,(double) i);

			if (d == 3) {
				GradChiRef_tmp[0][0] *=  2.0;
				GradChiRef_tmp[1][0] *=  2.0;
				GradChiRef_tmp[1][1] *=  2.0;
				GradChiRef_tmp[1][2] *=  2.0;
				GradChiRef_tmp[2][0] *=  2.0/3.0*sqrt(6.0)*a_n            ;
				GradChiRef_tmp[2][1] *=  1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)  ;
				GradChiRef_tmp[2][2] *= -1.0/6.0*sqrt(6.0)*(3.0*b_n+1.0)*i;
				GradChiRef_tmp[2][3] *=  1.0/2.0*sqrt(6.0)                ;
				GradChiRef_tmp[2][4] *= -1.0/2.0*sqrt(6.0)*(i+j)          ;
				if (i > 0) {
					GradChiRef_tmp[2][0] *= pow(1.0-b_n,i-1.0);
					GradChiRef_tmp[2][2] *= pow(1.0-b_n,i-1.0);
				} else {
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][2] *= 0.0;
				}
				GradChiRef_tmp[2][1] *= pow(1.0-b_n,(double) i);
				GradChiRef_tmp[2][3] *= pow(1.0-b_n,(double) i);
				GradChiRef_tmp[2][4] *= pow(1.0-b_n,(double) i);

				if (i+j > 0) {
					GradChiRef_tmp[0][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][1] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[1][2] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][0] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][1] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][2] *= pow(1.0-c_n,i+j-1.0);
					GradChiRef_tmp[2][4] *= pow(1.0-c_n,i+j-1.0);
				} else {
					GradChiRef_tmp[0][0] *= 0.0;
					GradChiRef_tmp[1][0] *= 0.0;
					GradChiRef_tmp[1][1] *= 0.0;
					GradChiRef_tmp[1][2] *= 0.0;
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][1] *= 0.0;
					GradChiRef_tmp[2][2] *= 0.0;
					GradChiRef_tmp[2][4] *= 0.0;
				}
				GradChiRef_tmp[2][3] *= pow(1.0-c_n,(double) (i+j));
			}

			// Sum contributions from all partial derivatives
			for (dim = 0; dim < d; dim++)
				for (pder = 0; pder < 5; pder++)
					GradChiRef_rst[dim][IndGrad] += GradChiRef_tmp[dim][pder];

			// Add scaling constant
			if (d == 2) {
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= 2.0/pow(3.0,0.25);
			} else {
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= 4.0/pow(2.0,0.25);
			}
		}
		Indbf++;
	}}}

	// Transpose GradChiRef_rst
	for (dim = 0; dim < d; dim++)
		mkl_dimatcopy('R','T',Nbf,Nn,1.,GradChiRef_rst[dim],Nn,Nbf);

//	array_print_d(Nn,Nbf,GradChiRef_rst[0],'R');

	free(a);
	free(b);
	free(c);

	*NbfOut = Nbf;
	return GradChiRef_rst;
}
