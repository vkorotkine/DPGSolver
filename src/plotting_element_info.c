#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <string.h>

//#include "database.h"
#include "parameters.h"
#include "functions.h"

//#include "petscsys.h"

/*
 *	Purpose:
 *		Return nodes and connectivity information required for plotting with paraview using the VTK format.
 *
 *	Comments:
 *		rst is stored in memory as r, s, then t.
 *		The ordering of 'Con' is determined by the by the cell type ordering in VTK for the following elements:
 *			QUAD, HEX, WEDGE, PYR.
 *
 *	Notation:
 *
 *	References:
 *		VTK User's Guide, File Formats for VTK Version 4.2 (Figure 2. Linear cell types found in VTK)
 */

void plotting_element_info(double **rst, unsigned int **connect, unsigned int **types, unsigned int *Nn,
                           const unsigned int P, const unsigned int typeIn)
{
	unsigned int i, j, k, l, m, iMax, jMax, kMax, lMax, row,
	             d, NnOut, NEOut,
	             *connectOut, *typesOut;
	double *rstOut;

	if (P == 0)
		printf("Error: Input P must be greater than 0 for plotting nodes.\n"), exit(1);

	if (typeIn == LINE || typeIn == QUAD || typeIn == HEX) {
		unsigned int dim,
		             N, u1,
		             Indr, Indc,
		             N2,
		             nLINE[2], nQUAD[4], nHEX[8];
		int          sd, sP, sN;
		double *r;

		// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
		d = 0;

		if      (typeIn == LINE) d = 1;
		else if (typeIn == QUAD) d = 2;
		else if (typeIn == HEX)  d = 3;

		N = P+1;
		NnOut = pow(N,d);
		NEOut = pow(P,d);

		u1 = 1;
		sd = d;
		sP = P;
		sN = N;

		r      = malloc(N       * sizeof *r);      // free
		rstOut = malloc(NnOut*d * sizeof *rstOut); // keep (requires external free)

		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut);   // keep (requires external free)

		if (P == 0) {
			r[0] = 0;
		} else {
			for (i = 0; i < N; i++)
				r[i] = -1.0 + (2.0*i)/P;
		}

		row = 0;
		for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
			for (dim = 0; dim < d; dim++) {
				if (dim == 0) Indr = i;
				if (dim == 1) Indr = j;
				if (dim == 2) Indr = k;
				rstOut[dim*NnOut+row] = r[Indr];
			}
			row++;
		}}}
		free(r);

		N2 = pow(N,2);

		Indc = 0;
		for (k = 0, kMax = (unsigned int) max(sP*min(sd-2,1),1); k < kMax; k++) {
			for (j = 0, jMax = max(P*min(d-1,u1),u1); j < jMax; j++) {
				for (i = 0; i < P; i++) {
					nLINE[0] = i;
					nLINE[1] = i+1;

					for (l = 0; l < 2; l++)             nQUAD[l] = nLINE[l]+N*j;
					for (l = 2, m = 1; l < 4; l++, m--) nQUAD[l] = nLINE[m] + N*(j+1);

					for (l = 0; l < 4; l++)             nHEX[l] = nQUAD[l] + N2*k;
					for (l = 4, m = 0; l < 8; l++, m++) nHEX[l] = nQUAD[m] + N2*(k+1);

					for (l = 0, lMax = pow(2,d); l < lMax; l++)
						connectOut[Indc*8+l] = nHEX[l];
					Indc++;
				}
			}
		}

		if (typeIn == LINE) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 3;
		} else if (typeIn == QUAD) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 9;
		} else if (typeIn == HEX) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 12;
		}

		*rst     = rstOut;
		*connect = connectOut;
		*types   = typesOut;
		*Nn      = NnOut;
	} else if (typeIn == TRI) {
		d = 2;
		NnOut = 1.0/2.0*(P+1)*(P+2);
		NEOut = 0;
		for (i = 0; i < P; i++)
			NEOut += 2*i+1;

		unsigned int iStart;
		double di, dj,
		       rst_c[3*d], BCoords[NnOut*3];

		rstOut     = malloc(NnOut*d * sizeof *rstOut); // keep (requires external free)
		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut); // keep (requires external free)

		// Determine barycentric coordinates of equally spaced nodes
		row = 0;
		for (i = 0; i <= P; i++) {
		for (j = 0, jMax = P-i; j <= jMax; j++) {
			di = i;
			dj = j;

			BCoords[row*3+2] = di/P;
			BCoords[row*3+1] = dj/P;
			BCoords[row*3+0] = 1.0 - (BCoords[row*3+1]+BCoords[row*3+2]);

			row++;
		}}

		// TRI corner nodes
		rst_c[0*2+0] = -1.0; rst_c[0*2+1] = -1.0/sqrt(3.0);
		rst_c[1*2+0] =  1.0; rst_c[1*2+1] = -1.0/sqrt(3.0);
		rst_c[2*2+0] =  0.0; rst_c[2*2+1] =  2.0/sqrt(3.0);

		rstOut = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NnOut,2,3,1.0,BCoords,rst_c);
		// keep (requires external free)

//array_print_d(NnOut,3,BCoords,'R');
//array_print_d(NnOut,d,rstOut,'R');

		row = 0;
		// Right-side up TRIs
		for (j = P; j; j--) {
			iStart = 0;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+1;
				connectOut[row*8+2] = i+1+j;
				row++;
			}
		}

		// Upside down TRIs
		for (j = P; j; j--) {
			iStart = 1;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j-1; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+j;
				connectOut[row*8+2] = i+j+1;
				row++;
			}
		}

//array_print_ui(NEOut,8,connectOut,'R');

		for (i = 0; i < NEOut; i++)
			typesOut[i] = 5;

		*rst     = rstOut;
		*connect = connectOut;
		*types   = typesOut;
		*Nn      = NnOut;

	} else if (typeIn == TET) {
		d = 3;
		NnOut = 1.0/6.0*(P+1)*(P+2)*(P+3);
		NEOut = 0;
		for (i = 1; i <= P; i++) {
		for (j = 1; j <= i; j++) {
			NEOut += j; // TETs
			if (i != P)
				NEOut += 2*j; // PYRs
		}}

		rstOut     = malloc(NnOut*d * sizeof *rstOut); // keep (requires external free)
		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut); // keep (requires external free)




		l = 0;
		for (i = P; i >= 1; i--) {
			if (i != P) {
				for (j = 1; j <= i; j++) {
				for (k = 0; k < 2*j; k++) {
					typesOut[l] = 14;
					l++;
				}}
			}
			for (j = 1; j <= i; j++) {
			for (k = 0; k < j; k++) {
				typesOut[l] = 10;
				l++;
			}}
		}

printf("%d\n",NEOut);
array_print_ui(1,NEOut,typesOut,'R');

		*rst     = rstOut;
		*connect = connectOut;
		*types   = typesOut;
		*Nn      = NnOut;
	} else if (typeIn == WEDGE) {
		printf("Add in support for WEDGE (plotting_element_info).\n"), exit(1);
	} else if (typeIn == PYR) {
		printf("Add in support for PYR (plotting_element_info).\n"), exit(1);
	}

//array_print_d(NnOut,d,*rst,'C');
//array_print_ui(NEOut,8,*connect,'R');
//array_print_ui(NEOut,1,*types,'R');

}
