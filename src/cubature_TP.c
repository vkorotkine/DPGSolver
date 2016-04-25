#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

#include "petscsys.h"
#include "mkl.h"

/*
 *	Purpose:
 *		Return nodes, weights and symmetries for tensor-product cubature depending on the nodetype.
 *
 *	Comments:
 *		Note that GL nodes are NOT the optimal integration nodes for d > 1. This was discussed in Hughes' FEM book and
 *		the pyFR code has some nodal sets available for low orders. Possibly implement this later (ToBeDeleted).
 *		The order of rst, w is important for minimizing memory stride while computing the length 2 Discrete Fourier
 *		Transform (ToBeModified).
 *		Ordering convention:
 *			GL/GLL: 2-blocks of -/+ node with abs(node) going from 1 to 0 followed by optional 0 node.
 *		rst is stored in memory as r, s, then t. This was motivated by the desire to be consistent with how XYZ
 *		coordinates are stored in each element, with this ordering being required so that vectorized operations can be
 *		performed on XYZ arrays.
 *
 *	Notation:
 *		rst   : Nodes array of dimension Nn*d (column-major storage)
 *		w     : Weights array of dimension Nn*1
 *		symms : Symmetries array
 *		        This always returns 1d symmetries in the current implementation (ToBeModified)
 *		Nn    : (N)umber of (n)odes
 *		Ns    : (N)umber of (s)ymmetries
 *
 *	References:
 *		GL  : http://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature/
 *		      content/GaussLegendre.m
 *		GLL : http://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights/content/
 *		      lglnodes.m
 */

void cubature_TP(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                 const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	// Standard datatypes
	unsigned int i, j, k, iMax, jMax, kMax, dim, u1,
	             N, rInd, row, Nrows, NsOut,
	             *symmsOut, *Indices;
	int          sd, sN;
	double       norm, swapd,
	             *r, *rold, *rdiff, *wOut, *r_d, *wOut_d, *r_std, *wOut_std,
	             *V, *a, *CM, *eigs;

	// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
	Nrows = 0;

	N = P+1;

	u1 = 1;
	sd = d;
	sN = N;

	r    = malloc(N * sizeof *r); // free
	wOut = malloc(N * sizeof *wOut); // free

	// Note: GLL must be first as "GL" is in "GLL"
	if (strstr(NodeType,"GLL") != NULL) {
		if (P == 0)
			printf("Error: Cannot use GLL nodes of order P0.\n"), exit(1);

		rold  = malloc(N * sizeof *rold); // free
		rdiff = malloc(N * sizeof *rdiff); // free

		// Use the Chebyshve-Guass-Lobatto nodes as the first guess
		for (i = 0; i < N; i++)
			r[i] = -cos(PI*(i)/P);
// array_print_d(1,N,r,'R');

		// Legendre Vandermonde Matrix
		V = malloc(N*N * sizeof *V); // free
		for (i = 0, iMax = N*N; i < iMax; i++)
			V[i] = 0.;

		/* Compute P_(N) using the recursion relation. Compute its first and second derivatives and update r using the
		Newton-Raphson method */
		for (i = 0; i < N; i++) rold[i] = 2.;
		for (i = 0; i < N; i++) rdiff[i] = r[i]-rold[i];
		norm = array_norm_d(N,rdiff,"Inf");

		while (norm > EPS) {
			for (i = 0; i < N; i++)
				rold[i] = r[i];
			for (j = 0; j < N; j++) {
				V[0*N+j] = 1.;
				V[1*N+j] = r[j];
			}

		for (i = 1; i < P; i++) {
		for (j = 0; j < N; j++) {
			V[(i+1)*N+j] = ((2*i+1)*r[j]*V[i*N+j] - i*V[(i-1)*N+j])/(i+1);
		}}

		for (j = 0; j < N; j++)
			r[j] = rold[j] - (r[j]*V[P*N+j]-V[(P-1)*N+j])/(N*V[P*N+j]);

		for (i = 0; i < N; i++)
			rdiff[i] = r[i]-rold[i];
		norm = array_norm_d(N,rdiff,"Inf");
		}

		for (j = 0; j < N; j++)
			wOut[j] = 2./(P*N*pow(V[P*N+j],2));

		free(rold);
		free(rdiff);
		free(V);

// array_print_d(1,N,r,'R');
// array_print_d(1,N,wOut,'R');
	} else if (strstr(NodeType,"GL") != NULL) {
		// Build the companion matrix CM
		/* CM is defined such that det(rI-CM)=P_n(r), with P_n(r) being the Legendre poynomial under consideration.
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

		// Determine the abscissas (r) and weights (w)
		/* Because det(rI-CM) = P_n(r), the abscissas are the roots of the characteristic polynomial, the eigenvalues of
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
			r[j]    = eigs[j];
			wOut[j] = 2*pow(CM[0*i+j],2);
		}
		free(CM);
		free(eigs);

// array_print_d(1,N,r,'R');
// array_print_d(1,N,wOut,'R');
	}

	// Re-arrange r and w for GL/GLL nodes
	r_std    = malloc(N * sizeof *r_std); // free
	wOut_std = malloc(N * sizeof *wOut_std); // free

	for (i = 0; i < N; i++) {
		r_std[i]    = r[i];
		wOut_std[i] = wOut[i];
	}

	if (N % 2 == 1) {
		k = (unsigned int) floor(N/2);

		r[N-1]    = r_std[k];
		wOut[N-1] = wOut_std[k];

		j = k;
		for (i = 0, iMax = N-1; i < iMax; ) {
			r[i]    = r_std[k-j];
			wOut[i] = wOut_std[k-j];
			i++;

			r[i]    = r_std[k+j];
			wOut[i] = wOut_std[k+j];
			i++;

			j--;
		}
	} else {
		k = N/2;
		j = k;
		for (i = 0, iMax = N; i < iMax; ) {
			r[i]    = r_std[k-j];
			wOut[i] = wOut_std[k-j];
			i++;

			r[i]    = r_std[k+j-1];
			wOut[i] = wOut_std[k+j-1];
			i++;

			j--;
		}
	}

	free(r_std);
	free(wOut_std);


	r_d    = malloc(pow(N,d)*d * sizeof *r_d);    // keep (requires external free)
	wOut_d = malloc(pow(N,d)   * sizeof *wOut_d); // free/keep (Conditional return_w)

	row = 0; Nrows = pow(N,d);
	for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
	for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
	for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
		wOut_d[row] = wOut[i];
		if (d == 2) wOut_d[row] *= wOut[j];
		if (d == 3) wOut_d[row] *= wOut[j]*wOut[k];
		for (dim = 0; dim < d; dim++) {
			if (dim == 0) rInd = i;
			if (dim == 1) rInd = j;
			if (dim == 2) rInd = k;
			r_d[dim*Nrows+row] = r[rInd];
		}
		row++;
	}}}
	free(r);
	free(wOut);

	// Compute symmetries
	NsOut = (unsigned int) ceil(N/2.0);
	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	if (N % 2 == 1) {
		for (i = 0; i < NsOut-1; i++)
			symmsOut[i] = 2;
		symmsOut[NsOut-1] = 1;
	} else {
		for (i = 0; i < NsOut; i++)
			symmsOut[i] = 2;
	}

	*rst   = r_d;
	*symms = symmsOut;

	*Nn = pow(N,d);
	*Ns = NsOut;

	if (return_w != 0) *w = wOut_d;
	else               *w = NULL, free(wOut_d);

// array_print_d(pow(N,d),d,r_d,'R');
// array_print_d(pow(N,d),1,wOut_d,'R');
}
