#include <stdlib.h>
#include <stdio.h>

#include "parameters.h"
#include "functions.h"

#include "petscsys.h"
#include "mkl.h"

/*
 *	Purpose:
 *		Return nodes and weights for tensor-product cubature depending on the nodetype.
 *
 *	Comments:
 *		Pointers are only returned to desired variables (as indicated in ToReturn); other variables are freed.
 *			Note: xir and Nn are always returned
 *		Possibly modify ordering in connect based on vtk formatting convention. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 *		GL  : http://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature/
 *		      content/GaussLegendre.m
 *		GLL : http://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights/content/
 *		      lglnodes.m
 */

void cubature_TP(double **xir, double **W, unsigned int **Con, unsigned int *Nn, const unsigned int *ToReturn,
                 const unsigned int P, const unsigned int d,const char *NodeType)
{

	// Standard datatypes
	unsigned int i, j, k, iMax, jMax, kMax, dim, l, lMax, m, u1,
	             N, N2, xInd, row, nLINE[2], nQUAD[4], nHEX[8], IndC,
	             *Indices, *connect;
	int          sd, sP, sN;
	double       norm, swapd,
	             *x, *xold, *xdiff, *w, *x_d, *w_d,
	             *V, *a, *CM, *eigs;

	N = P+1;

	u1 = 1;
	sd = d;
	sP = P;
	sN = N;

	x = malloc(N * sizeof *x); // free
	w = malloc(N * sizeof *w); // free

	// Note: GLL must be first as "GL" is in "GLL"
	if (strstr(NodeType,"GLL") != NULL) {
		if (P == 0)
			printf("Error: Cannot use GLL nodes of order P0.\n"), exit(1);

		xold  = malloc(N * sizeof *xold); // free
		xdiff = malloc(N * sizeof *xdiff); // free

		// Use the Chebyshve-Guass-Lobatto nodes as the first guess
		for (i = 0; i < N; i++)
			x[i] = -cos(PI*(i)/P);
// array_print_d(1,N,x,'R');

		// Legendre Vandermonde Matrix
		V = malloc(N*N * sizeof *V); // free
		for (i = 0, iMax = N*N; i < iMax; i++)
			V[i] = 0.;

		/* Compute P_(N) using the recursion relation. Compute its first and second derivatives and update x using the
		Newton-Raphson method */
		for (i = 0; i < N; i++) xold[i] = 2.;
		for (i = 0; i < N; i++) xdiff[i] = x[i]-xold[i];
		norm = array_norm_d(N,xdiff,"Inf");

		while (norm > EPS) {
			for (i = 0; i < N; i++)
				xold[i] = x[i];
			for (j = 0; j < N; j++) {
				V[0*N+j] = 1.;
				V[1*N+j] = x[j];
			}

		for (i = 1; i < P; i++) {
		for (j = 0; j < N; j++) {
			V[(i+1)*N+j] = ((2*i+1)*x[j]*V[i*N+j] - i*V[(i-1)*N+j])/(i+1);
		}}

		for (j = 0; j < N; j++)
			x[j] = xold[j] - (x[j]*V[P*N+j]-V[(P-1)*N+j])/(N*V[P*N+j]);

		for (i = 0; i < N; i++)
			xdiff[i] = x[i]-xold[i];
		norm = array_norm_d(N,xdiff,"Inf");
		}

		for (j = 0; j < N; j++)
			w[j] = 2./(P*N*pow(V[P*N+j],2));

		free(xold);
		free(xdiff);
		free(V);

// array_print_d(1,N,x,'R');
// array_print_d(1,N,w,'R');
	} else if (strstr(NodeType,"GL") != NULL) {
		// Build the companion matrix CM
		/* CM is defined such that det(xI-CM)=P_n(x), with P_n(x) being the Legendre poynomial under consideration.
		 * Moreover, CM is constructed in such a way so as to be symmetrical.
		 */
		a = malloc(P * sizeof *a); // free
		for (i = 1; i < N; i++)
			a[i-1] = (1.*i)/sqrt(4.*pow(i,2)-1);

		CM = malloc(N*N * sizeof *CM); // free
		for (i = 0, iMax = N*N; i < iMax; i++)
			CM[i] = 0.;

		for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j+1) CM[i*N+j] = a[j];
			if (j == i+1) CM[i*N+j] = a[i];
		}}
		free(a);

		// Determine the abscissas (x) and weights (w)
		/* Because det(xI-CM) = P_n(x), the abscissas are the roots of the characteristic polynomial, the eigenvalues of
		 * CM. The weights can then be derived from the corresponding eigenvectors.
		 */

		eigs    = malloc(N *sizeof *eigs); // free
		Indices = malloc(N *sizeof *Indices); // free
		for (i = 0; i < N; i++)
			Indices[i] = i;

		if (LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',(MKL_INT) N,CM,(MKL_INT) N,eigs) > 0)
			printf("Error: mkl LAPACKE_sysev failed to compute eigenvalues.\n"), exit(1);

// array_print_d(1,N,eigs,'R');
// array_print_d(N,N,CM,'R');

		array_sort_d(1,N,eigs,Indices,'R','N');

		for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			swapd = CM[i*N+j];
			CM[i*N+j]          = CM[Indices[i]*N+j];
			CM[Indices[i]*N+j] = swapd;
		}}
		free(Indices);

		for (j = 0; j < N; j++) {
			x[j] = eigs[j];
			w[j] = 2*pow(CM[0*i+j],2);
		}
		free(CM);
		free(eigs);

// array_print_d(1,N,x,'R');
// array_print_d(1,N,w,'R');
	}

	x_d     = malloc(pow(N,d)*d        * sizeof *x_d); // keep
	w_d     = malloc(pow(N,d)*d        * sizeof *w_d); // free/keep (Conditional ToReturn)
	connect = malloc(pow(P,d)*pow(2,d) * sizeof *connect); // free/keep (Conditional ToReturn)

	if (strstr(NodeType,"GLL") != NULL || strstr(NodeType,"GL") != NULL) {
		row = 0;
		for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
			w_d[row] = w[i];
			if (d == 2) w_d[row] *= w[j];
			if (d == 3) w_d[row] *= w[j]*w[k];
			for (dim = 0; dim < d; dim++) {
				if (dim == 0) xInd = i;
				if (dim == 1) xInd = j;
				if (dim == 2) xInd = k;
				x_d[row*d+dim] = x[xInd];
			}
			row++;
		}}}
		free(x);
		free(w);

		*xir = x_d;
		*Nn  = pow(N,d);

		if (ToReturn[1] != 0) *W = w_d;
		else                  *W = NULL, free(w_d);

		if (ToReturn[2] != 0) {
			printf("Error: Connectivity can only be returned for the \"ES\" nodetype.\n"), exit(1);
		} else {
			*Con = NULL;
			free(connect);
		}

// array_print_d(pow(N,d),d,x_d,'R');
// array_print_d(pow(N,d),1,w_d,'R');

	return;
	}

	if (strstr(NodeType,"ES") != NULL) {
		free(w);
		if (P == 0) {
			x[0] = 0;
		} else {
			for (i = 0; i < N; i++)
				x[i] = -1.+2./P*i;
		}

		row = 0;
		for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
			for (dim = 0; dim < d; dim++) {
				if (dim == 0) xInd = i;
				if (dim == 1) xInd = j;
				if (dim == 2) xInd = k;
				x_d[row*d+dim] = x[xInd];
			}
			row++;
		}}}
		free(x);

		N2 = pow(N,2);
		for (i = 0, IndC = 0; i < P; i++) {
			nLINE[0] = i;
			nLINE[1] = i+1;
			for (j = 0, jMax = max(P*min(d-1,u1),u1); j < jMax; j++) {
				for (l = 0; l < 2; l++)             nQUAD[l] = nLINE[l]+N*j;
				for (l = 2, m = 1; l < 4; l++, m--) nQUAD[l] = nLINE[m] + N*(j+1);
				for (k = 0, kMax = (unsigned int) max(sP*min(sd-2,1),1); k < kMax; k++) {
					for (l = 0; l < 4; l++)             nHEX[l] = nQUAD[l] + N2*k;
					for (l = 4, m = 0; l < 8; l++, m++) nHEX[l] = nQUAD[m] + N2*(k+1);
					for (l = 0, lMax = pow(2,d); l < lMax; l++)
						connect[IndC*lMax+l] = nHEX[l];
						IndC++;
				}
			}
		}

// array_print_i(pow(P,d),pow(2,d),connect,'R');

		*xir = x_d;
		*Nn  = pow(N,d);
		if (ToReturn[1] != 0) {
			printf("Error: There are no weights associated with the \"ES\" nodetype.\n"), exit(1);
		} else {
			*W = NULL;
			free(w_d);
		}

		if (ToReturn[2] != 0) *Con = connect;
		else                  *Con = NULL, free(connect);
	}

}
