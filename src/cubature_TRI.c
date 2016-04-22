#include <stdlib.h>
#include <stdio.h>
//#include <math.h>

//#include "parameters.h"
//#include "functions.h"

//#include "petscsys.h"
//#include "mkl.h"

/*
 *	Purpose:
 *		Return nodes, weights and symmetries for triangular cubature depending on the nodetype.
 *
 *	Comments:
 *		Check that AO nodes correspond to those in pyfr after writing the converter script. (ToBeDeleted)
 *		The WS and WV nodes were determined based off of those from the pyfr code (pyfr/quadrules/tri) after being
 *		transfered to the equilateral reference TRI used in this code.
 *		The order of r, w is important for minimizing memory stride while computing the length 3 Discrete Fourier
 *		Transform (ToBeModified).
 *		Ordering convention:
 *			3-blocks of symmetric nodes going from farthest from center towards the center, followed by 1-block of
 *			center node if present.
 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).

 *	Notation:
 *		rst   : Nodes array of dimension Nn*d (column-major storage)
 *		w     : Weights array of dimension Nn*1
 *		symms : Symmetries array
 *		        This always returns 1d symmetries in the current implementation (ToBeModified)
 *		Nn    : (N)umber of (n)odes
 *		Ns    : (N)umber of (s)ymmetries
 *
 *	References:
 *		pyfr code : http://www.pyfr.org
 *
 *		AO : Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods
 *		WS : ToBeModified
 *		WV : ToBeModified
 */

void cubature_TRI(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d,const char *NodeType)
{
	// Standard datatypes
	unsigned int *symmsOut;
	double       *rstOut, *wOut;

	// Silence compiler warnings
	Nrows = 0;

/*
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
*/

	*rst   = rOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w != 0) *w = wOut;
	else               *w = NULL, free(wOut);
}
