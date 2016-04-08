#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Return the "matrix" GradChiRef_xir representing the gradients of the orthonormal basis functions evaluated at
 *		the provided quadrature nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrix" is returned as d pointers to 1D arrays.
 *
 *	Notation:
 *
 *	Example code:
 *
 *		int P_tmp, Nn_tmp, dim_tmp,
 *		    *dummyi_tmp, ToReturn_tmp[4];
 *		double *xir_tmp, *dummyd_tmp, **GradChiRef_tmp;
 *
 *		P_tmp = 2;
 *		dim_tmp = 2;
 *
 *		ToReturn_tmp[0] = 1; ToReturn_tmp[1] = 0; ToReturn_tmp[2] = 0; ToReturn_tmp[3] = 1;
 *		cubature_TP(&xir_tmp,&dummyd_tmp,&dummyi_tmp,&Nn_tmp,ToReturn_tmp,P_tmp,dim_tmp,"GL");
 *		GradChiRef_tmp = grad_basis_TP(P_tmp,xir_tmp,Nn_tmp,dim_tmp);
 *
 *		array_print_d(Nn_tmp,pow(P_tmp+1,dim_tmp),GradChiRef_tmp[0],'R');
 *
 *		free(xir_tmp);
 *		array_free2_d(dim_tmp,GradChiRef_tmp);
 *
 *	References:
 *
 */

double **grad_basis_TP(const int P, const double *xir, const int Nn, const int d)
{
	int    i, j, k, iMax, jMax, kMax, dim, Indbf, Indn,
	       N, Nbf;
	double **GradChiRef_xir;

	N = P+1;
	Nbf = pow(N,d);

	GradChiRef_xir = malloc(d * sizeof *GradChiRef_xir); // keep
	for (dim = 0; dim < d; dim++)
		GradChiRef_xir[dim] = malloc(Nn*Nbf * sizeof **GradChiRef_xir); // keep

	Indbf = 0;
	if (d == 1) {
		for (i = 0, iMax = min(max((d)  *N,1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_xir[0][Indbf*Nn+Indn] = grad_jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i);
			}
			Indbf++;
		}
	} else if (d == 2) {
		for (j = 0, jMax = min(max((d-1)*N,1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d)  *N,1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_xir[0][Indbf*Nn+Indn] = grad_jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i)*
				                                        jacobiP(xir[Indn*d+1],0.0,0.0,(double) j);
				GradChiRef_xir[1][Indbf*Nn+Indn] =      jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i)*
				                                   grad_jacobiP(xir[Indn*d+1],0.0,0.0,(double) j);
			}
			Indbf++;
		}}
	} else if (d == 3) {
		for (k = 0, kMax = min(max((d-2)*N,1),N); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d)  *N,1),N); i < iMax; i++) {
			for (Indn = 0; Indn < Nn; Indn++) {
				GradChiRef_xir[0][Indbf*Nn+Indn] = grad_jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i)*
				                                        jacobiP(xir[Indn*d+1],0.0,0.0,(double) j)*
				                                        jacobiP(xir[Indn*d+2],0.0,0.0,(double) k);
				GradChiRef_xir[1][Indbf*Nn+Indn] =      jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i)*
				                                   grad_jacobiP(xir[Indn*d+1],0.0,0.0,(double) j)*
				                                        jacobiP(xir[Indn*d+2],0.0,0.0,(double) k);
				GradChiRef_xir[2][Indbf*Nn+Indn] =      jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i)*
				                                        jacobiP(xir[Indn*d+1],0.0,0.0,(double) j)*
				                                   grad_jacobiP(xir[Indn*d+2],0.0,0.0,(double) k);
			}
			Indbf++;
		}}}
	}

	// Transpose GradChiRef_xir
	for (dim = 0; dim < d; dim++)
		mkl_dimatcopy('R','T',Nbf,Nn,1.,GradChiRef_xir[dim],Nn,Nbf);

// array_print_d(Nn,Nbf,GradChiRef_xir[0],'R');

	return GradChiRef_xir;
}
