#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

#include "mkl.h"

/*
 *	Purpose:
 *		Return the "matrix" ChiRef_xir representing the orthonormal basis functions evaluated at the provided quadrature
 *		nodes for polynomial order P.
 *
 *	Comments:
 *		The "matrix" is returned as a 1D array.
 *
 *	Notation:
 *
 *	Example code:
 *
 *		int P_tmp, Nn_tmp, dim_tmp,
 *		    *dummyi_tmp, ToReturn_tmp[4];
 *		double *xir_tmp, *dummyd_tmp, *ChiRef_tmp;
 *
 *		P_tmp = 2;
 *		dim_tmp = 2;
 *
 *		ToReturn_tmp[0] = 1; ToReturn_tmp[1] = 0; ToReturn_tmp[2] = 0; ToReturn_tmp[3] = 1;
 *		cubature_TP(&xir_tmp,&dummyd_tmp,&dummyi_tmp,&Nn_tmp,ToReturn_tmp,P_tmp,dim_tmp,"GL");
 *		ChiRef_tmp = basis_TP(P_tmp,xir_tmp,Nn_tmp,dim_tmp);
 *
 *		array_print_d(Nn_tmp,pow(P_tmp+1,dim_tmp),ChiRef_tmp);
 *
 *		free(xir_tmp);
 *		free(ChiRef_tmp);
 *
 *	References:
 *
 */

double *basis_TP(const int P, const double *xir, const int Nn, const int d)
{
	int    i, j, k, iMax, jMax, kMax, Indbf, Indn,
	       N, Nbf;
	double *ChiRef_xir;

    N = P+1;
	Nbf = pow(N,d);

	ChiRef_xir = malloc(Nn*Nbf * sizeof *ChiRef_xir); // keep

	Indbf = 0;
	for (k = 0, kMax = min(max((d-2)*N,1),N); k < kMax; k++) {
	for (j = 0, jMax = min(max((d-1)*N,1),N); j < jMax; j++) {
	for (i = 0, iMax = min(max((d)  *N,1),N); i < iMax; i++) {
		for (Indn = 0; Indn < Nn; Indn++) {
			           ChiRef_xir[Indbf*Nn+Indn] =  jacobiP(xir[Indn*d]  ,0.0,0.0,(double) i);
			if (d > 1) ChiRef_xir[Indbf*Nn+Indn] *= jacobiP(xir[Indn*d+1],0.0,0.0,(double) j);
			if (d > 2) ChiRef_xir[Indbf*Nn+Indn] *= jacobiP(xir[Indn*d+2],0.0,0.0,(double) k);
		}
		Indbf++;
	}}}

	// Transpose ChiRef_xir
	mkl_dimatcopy('R','T',Nbf,Nn,1.,ChiRef_xir,Nn,Nbf);

// array_print_d(Nn,Nbf,ChiRef_xir);

	return ChiRef_xir;
}
