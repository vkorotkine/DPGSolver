// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

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
 */

double **grad_basis_PYR(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                        const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, pder, Indbf, Indn, IndGrad, Nbf, mu_ij;
	double       **GradChiRef_rst, *a, *b, *c, a_n, b_n, c_n, jPa, jPb, jPc, djPa, djPb, djPc;
	double       GradChiRef_tmp[3][5];

	if (d != 3)
		printf("Error: grad_basis_PYR only supports d = 3.\n"), exit(1);

	// Convert from rst to abc coordinates
	a = malloc(Nn * sizeof *a); // free
	b = malloc(Nn * sizeof *b); // free
	c = malloc(Nn * sizeof *c); // free

	rst_to_abc_PYR(Nn,d,rst,a,b,c);

	Nbf = (unsigned int) (1.0/6.0*((P+1)*(P+2)*(2*P+3)));

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		 GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P; i <= iMax; i++) {
	for (j = 0, jMax = P; j <= jMax; j++) {
		mu_ij = max(i,j);
		for (k = 0, kMax = P-mu_ij; k <= kMax; k++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				IndGrad = Indbf*Nn+Indn;
				a_n = a[Indn];
				b_n = b[Indn];
				c_n = c[Indn];

				jPa = jacobiP(a_n,0.0,            0.0,(double) i);
				jPb = jacobiP(b_n,0.0,            0.0,(double) j);
				jPc = jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

				djPa = grad_jacobiP(a_n,0.0,            0.0,(double) i);
				djPb = grad_jacobiP(b_n,0.0,            0.0,(double) j);
				djPc = grad_jacobiP(c_n,2.0*(mu_ij+1.0),0.0,(double) k);

				for (dim = 0; dim < d; dim++)
					for (pder = 0; pder < 4; pder++ )
						GradChiRef_tmp[dim][pder] = 0.0;

				// Obtain contributions from each partial derivative
				GradChiRef_tmp[0][0] = djPa* jPb* jPc;
				GradChiRef_tmp[1][1] =  jPa*djPb* jPc;
				GradChiRef_tmp[2][0] = djPa* jPb* jPc;
				GradChiRef_tmp[2][1] =  jPa*djPb* jPc;
				GradChiRef_tmp[2][2] =  jPa* jPb*djPc;
				GradChiRef_tmp[2][3] =  jPa* jPb* jPc;

				GradChiRef_tmp[0][0] *=  2.0          ;
				GradChiRef_tmp[1][1] *=  2.0          ;
				GradChiRef_tmp[2][0] *=  sqrt(2.0)*a_n;
				GradChiRef_tmp[2][1] *=  sqrt(2.0)*b_n;
				GradChiRef_tmp[2][2] *=  sqrt(2.0)    ;
				GradChiRef_tmp[2][3] *= -sqrt(2.0)*mu_ij    ;

				if (mu_ij > 0) {
					GradChiRef_tmp[0][0] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[1][1] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][0] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][1] *= pow(1.0-c_n,mu_ij-1.0);
					GradChiRef_tmp[2][3] *= pow(1.0-c_n,mu_ij-1.0);
				} else {
					GradChiRef_tmp[0][0] *= 0.0;
					GradChiRef_tmp[1][1] *= 0.0;
					GradChiRef_tmp[2][0] *= 0.0;
					GradChiRef_tmp[2][1] *= 0.0;
					GradChiRef_tmp[2][3] *= 0.0;
				}
				GradChiRef_tmp[2][2] *= pow(1.0-c_n,(double) mu_ij);

				// Sum contributions from all partial derivatives
				for (dim = 0; dim < d; dim++)
					for (pder = 0; pder < 4; pder++)
						GradChiRef_rst[dim][IndGrad] += GradChiRef_tmp[dim][pder];

				// Add scaling constant
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][IndGrad] *= pow(2.0,1.25);
			}
			Indbf++;
		}
	}}

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
