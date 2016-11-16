// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "solver_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Parameters.h"

/*
 *	Purpose:
 *		Provide solver-related functions.
 *
 *	Comments:
 *		If it is found that a significant portion of time is spent in this function after profiling, compare speed with
 *		the addition of switch statement with explicitly written out low-order options. (ToBeDeleted)
 *		TRIs are significantly more complicated to handle because the the somewhat arbitrary position of the nodes
 *		within the 3-symmetry orbits. Thus, the physical coordinates are used to swap indices for the last three cases,
 *		which ensures that the ordering will be correct as long as nodes are input in 3 and 1-symmetry blocks.
 *
 *	Notation:
 *
 *	References:
 */

void get_face_ordering(const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
                        const unsigned int Nn, const unsigned int Ns, const unsigned int *symms, const double *rst,
                        unsigned int *nOrd)
{
	/*
	 *	Purpose:
	 *		Return ordering of opposite FACE such that surface nodes match.
	 *
	 *	Comments:
	 *		FType is only used in 3D.
	 *		symms is only used for 3D TRI FACEs.
	 */

	unsigned int i, j, iMax, jMax, iInd;

	switch (d) {
	case 1:
		nOrd[0] = 0;
		break;
	case 2:
		switch (IndOrd) {
		default:
			for (i = 0; i < Nn; i++)
				nOrd[i] = i;
			break;
		case 1:
			// Add in switch (Nn) here and write out low order options (ToBeDeleted)

			if (Nn % 2 == 0) {
				for (i = 0, iMax = Nn; i < iMax; i++) {
					if (i % 2 == 0) nOrd[i] = i+1;
					else            nOrd[i] = i-1;
				}
			} else {
				for (i = 0, iMax = Nn-1; i < iMax; i++) {
					if (i % 2 == 0) nOrd[i] = i+1;
					else            nOrd[i] = i-1;
				}
				nOrd[iMax] = iMax;
			}
			break;
		}
		break;
	default: // default to 3D
		if (FType == QUAD) {
			unsigned int sqrtNn = sqrt(Nn), nOrdswap[Nn];
			// Add in switch (Nn) here and write out low order options (ToBeDeleted)

			switch(IndOrd) {
				default: // default case 0
				case 1:
				case 2:
				case 3:
					for (i = 0; i < Nn; i++)
						nOrd[i] = i;
					break;
				case 4:
				case 5:
				case 6:
				case 7:
					for (i = 0; i < sqrtNn; i++) {
						iInd = i*sqrtNn;
						for (j = 0; j < sqrtNn; j++)
							nOrd[iInd+j] = i+j*sqrtNn;
					}
					break;
			}

			// Swap 1D-blocks if necessary
			switch(IndOrd) {
				default:
					; // Do nothing
					break;
				case 2:
				case 3:
				case 6:
				case 7:
					for (i = 0; i < Nn; i++)
						nOrdswap[i] = nOrd[i];

					if (sqrtNn % 2 == 0) iMax = sqrtNn;
					else                 iMax = sqrtNn-1;

					for (i = 0; i < iMax; i++) {
						iInd = i*sqrtNn;
						for (j = 0, jMax = sqrtNn; j < jMax; j++) {
							if (i % 2 == 0) nOrd[iInd+j] = nOrdswap[(iInd+sqrtNn)+j];
							else            nOrd[iInd+j] = nOrdswap[(iInd-sqrtNn)+j];
						}
					}
					// Setting the last block of nOrd is redundant for sqrtNn odd as it is unchanged.
					break;
			}

			// Reverse entries of 1D-blocks if necessary
			switch(IndOrd) {
				default:
					return;
					break;
				case 1:
				case 3:
				case 5:
				case 7:
					for (i = 0; i < Nn; i++)
						nOrdswap[i] = nOrd[i];

					if (sqrtNn % 2 == 0) jMax = sqrtNn;
					else                 jMax = sqrtNn-1;

					for (i = 0, iMax = sqrtNn; i < iMax; i++) {
						iInd = i*iMax;
						for (j = 0; j < jMax; j++) {
							if (j % 2 == 0) nOrd[iInd+j] = nOrdswap[iInd+j+1];
							else            nOrd[iInd+j] = nOrdswap[iInd+j-1];
						}
					}
					break;
			}
		} else if (FType == TRI) {
			unsigned int j, k, kMax, iInd, subOrder[3], nOrdswap3[3], nOrdswap[Nn], Foundn[Nn], IndX[Nn];
			double       DY[Nn*Nn];
			// Add in switch (Nn) here and write out low order options (ToBeDeleted)

			for (i = 0; i < Nn; i++)
				nOrd[i] = i;

			// Swap entries if necessary
			switch(IndOrd) {
				default: // default cases 0, 1, 2
					; // Do nothing
					break;
				case 3:
				case 4:
				case 5:
					for (i = 0; i < Nn; i++) {
						iInd = i*Nn;
						for (j = 0; j < Nn; j++) {
							DY[iInd+j] = fabs(rst[Nn+i]-rst[Nn+j]);
						}
					}

					for (i = 0; i < Nn; i++)
						Foundn[i] = 0;

					for (i = 0; i < Nn; i++) {
						if (!Foundn[i]) {
							iInd = i*Nn;
							kMax = 0;
							for (j = 0; j < Nn; j++) {
								if (!Foundn[j] && i != j && DY[iInd+j] < 10*EPS)
									IndX[kMax++] = j;
							}
							for (k = 0; k < kMax; k++) {
								if (fabs(rst[i]+rst[IndX[k]]) < 1e3*EPS) {
									Foundn[i] = 1;
									Foundn[IndX[k]] = 1;
									nOrdswap[i] = IndX[k];
									nOrdswap[IndX[k]] = i;

									break;
								}
							}
							if (kMax == 0) {
								Foundn[i] = 1;
								nOrdswap[i] = i;
							}
						}
					}

					for (i = 0; i < Nn; i++) {
						if (Foundn[i] == 0)
							printf("Error: Did not find all nodes in get_face_ordering (TRI).\n"), exit(1);
					}

					for (i = 0; i < Nn; i++)
						nOrd[i] = nOrdswap[i];
					break;
			}

			// Rotate entries of 3-symmetry blocks if necessary
			switch(IndOrd) {
				default: // cases 0, 5
					; // No rotations needed
					return;
					break;
				case 1:
				case 3:
					subOrder[0] = 1; subOrder[1] = 2; subOrder[2] = 0;
					break;
				case 2:
				case 4:
					subOrder[0] = 2; subOrder[1] = 0; subOrder[2] = 1;
					break;
			}

			iInd = 0;
			for (i = 0; i < Ns; i++) {
				if (i) iInd += symms[i-1];
				jMax = symms[i];
				if (jMax == 3) {
					for (j = 0; j < jMax; j++)
						nOrdswap3[j] = nOrd[iInd+subOrder[j]];
					for (j = 0; j < jMax; j++)
						nOrd[iInd+j] = nOrdswap3[j];
				}
				// Setting the 1-symmetry orbit (if present) is redundant as the node position remains unchanged.
			}
		} else {
			printf("Error: Unsupported FType in 3D in get_face_ordering.\n"), exit(1);
		}
		break;
	}
}
