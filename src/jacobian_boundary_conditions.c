// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "jacobian_boundary_conditions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "S_DB.h"

#include "variable_functions.h"

/*
 *	Purpose:
 *		Compute the jacobian of boundary conditions for given input parameters: XYZ, WL, WOut, nL
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

/*
void jacobian_boundary_Riemann(const unsigned int Nn, const unsigned int Nel, double *XYZ, double *WL, double *WOut,
                               double *WB, double *nL, const unsigned int d)
{

}
*/

void jacobian_boundary_SlipWall(const unsigned int Nn, const unsigned int Nel, double *WL, double *dWdW, double *nL,
                                const unsigned int d)
{
	// Standard datatypes
	unsigned int i, NnTotal, IndE;
	double *rhoL, *rhouL, *rhovL, *rhowL, *EL, *rhoB, *rhouB, *rhovB, *rhowB, *EB, rhoVL;

	NnTotal = Nn*Nel;
	IndE = d+1;

	if (d == 3) {
	} else if (d == 2) {
	} else if (d == 1) {
	}
}
