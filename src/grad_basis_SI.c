#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Return the "matrices" ChiRef_rst representing the gradients of the orthonormal basis functions evaluated at the
 *		provided quadrature nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrices" are returned as a 1D array.
 *
 *	Notation:
 *
 *	References:
 */

double **grad_basis_SI(const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                       const unsigned int d)
{
	unsigned int i, j, k, iMax, jMax, kMax, dim, Indbf, Indn, IndGrad, Nbf;
	double       **GradChiRef_rst, *a, *b, *c;

	if (d < 2 || d > 3)
		printf("Error: grad_basis_SI only supports d = [2,3].\n"), exit(1);

	// Convert from rst to abc coordinates
	a = calloc(Nn , sizeof *a); // free
	b = calloc(Nn , sizeof *b); // free
	c = calloc(Nn , sizeof *c); // free

	rst_to_abc(Nn,d,rst,a,b,c);

	Nbf = factorial_ui(d+P)/(factorial_ui(d)*factorial_ui(P));

	GradChiRef_rst = malloc(d * sizeof *GradChiRef_rst); // keep (requires external free)
	for (dim = 0; dim < d; dim++)
		 GradChiRef_rst[dim] = calloc(Nn*Nbf , sizeof **GradChiRef_rst); // keep (requires external free)
//		 GradChiRef_rst[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_rst); // keep (requires external free)

	Indbf = 0;
	for (i = 0, iMax = P;   i <= iMax; i++) {
	for (j = 0, jMax = P-i; j <= jMax; j++) {
	for (k = 0, kMax = (d-2)*(P-i-j); k <= kMax; k++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			IndGrad = Indbf*Nn+Indn;
			a_n = a[Indn];
			b_n = b[Indn];
			c_n = c[Indn];

			grad_a_n = grad_a[Indn*3];
			grad_b_n = grad_b[Indn*3];
			grad_c_n = grad_c[Indn*3];

			for (dim = 0; dim < d; dim++)
				GradChiRef_rst[0][IndGrad] = 0.0;

			if (i > 0) {
				GradChiRef_rst[1][IndGrad] += i*pow(1.0-b_n,i-1.0)*grad_b_n[1];
				                             *jacobiP(a_n,0.0,      0.0,(double i))
				                             *jacobiP(b_n,2.0*i+1.0,0.0,(double j));
				//GradChiRef_rst[2][IndGrad]
			}

			GradChiRef_rst[0][IndGrad] += pow(1.0-b_n,i)
			                             *





			if (d == 2) {
				for (dim = 0; dim < d; dim++)
					GradChiRef_rst[dim][Indbf*Nn+Indn] *= 2.0/pow(3.0,0.25);
			} else {
				for (dim = 0; dim < 2; dim++)
					GradChiRef_rst[dim][Indbf*Nn+Indn] *= 2.0/pow(3.0,0.25)*pow(1.0-c[Indn],i+j);
				dim = 2;
				GradChiRef_rst[dim][Indbf*Nn+Indn] *= 2.0/pow(3.0,0.25);
			}



/*
			if (d == 2)
				ChiRef_rst[Indbf*Nn+Indn] = 2.0/pow(3.0,0.25)*pow(1.0-b[Indn],i);
			else
				ChiRef_rst[Indbf*Nn+Indn] = 4.0/pow(2.0,0.25)*pow(1.0-b[Indn],i)*pow(1.0-c[Indn],i+j);

			ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(a[Indn],0.0,0.0,(double) i);
			ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(b[Indn],2.0*i+1.0,0.0,(double) j);
			if (d == 3)
				ChiRef_rst[Indbf*Nn+Indn] *= jacobiP(c[Indn],2.0*(i+j+1.0),0.0,(double) k);
*/
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
