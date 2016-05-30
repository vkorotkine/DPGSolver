// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Provide simple solver-related functions:
 *			void get_facet_ordering(const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
 *			                        const unsigned int Nn, const unsigned int Ns, const unsigned int symms,
 *			                        unsigned int nOrd)
 *
 *	Comments:
 *		If it is found that a significant portion of time is spent in this function after profiling, compare speed with
 *		the addition of switch statement with explicitly written out low-order options. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

void get_facet_ordering(const unsigned int d, const unsigned int IndOrd, const unsigned int FType,
                        const unsigned int Nn, const unsigned int Ns, const unsigned int *symms, unsigned int *nOrd)
{
	/*
	 *	Purpose:
	 *		Return ordering of opposite FACET such that surface nodes match.
	 *
	 *	Comments:
	 *		FType is only used in 3D.
	 *		symms is only used for 3D TRI FACETs.
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
			unsigned int subOrder[3];
			// Add in switch (Nn) here and write out low order options (ToBeDeleted)

			switch(IndOrd) {
			default: // default case 0
				subOrder[0] = 0; subOrder[1] = 1; subOrder[2] = 2;
				break;
			case 1:
				subOrder[0] = 1; subOrder[1] = 2; subOrder[2] = 0;
				break;
			case 2:
				subOrder[0] = 2; subOrder[1] = 0; subOrder[2] = 1;
				break;
			case 3:
				subOrder[0] = 0; subOrder[1] = 2; subOrder[2] = 1;
				break;
			case 4:
				subOrder[0] = 2; subOrder[1] = 1; subOrder[2] = 0;
				break;
			case 5:
				subOrder[0] = 1; subOrder[1] = 0; subOrder[2] = 2;
				break;
			}

			iInd = 0;
			for (i = 0; i < Ns; i++) {
				if (i) iInd += symms[i-1];
				jMax = symms[i];
				if (jMax == 3) {
					for (j = 0; j < jMax; j++) {
						nOrd[iInd+j] = iInd+subOrder[j];
					}
				} else { // jMax == 1
					nOrd[iInd] = iInd;
				}
			}
		} else {
			printf("Error: Unsupported FType in 3D in get_facet_ordering.\n"), exit(1);
		}
		break;
	}


}
