// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_Advection_functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

/*
 *	Purpose:
 *		Provide Advection related solver functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

double *compute_b_Advection(unsigned int const Nn, double const *const XYZ)
{
	char const *const PDESpecifier = DB.PDESpecifier;
	unsigned int const d = DB.d;

	double *const b = malloc(Nn*d * sizeof *b); // keep

	if (strstr(PDESpecifier,"Default")) {
		if (strstr(PDESpecifier,"Unsteady")) {
			EXIT_UNSUPPORTED;
		} else if (strstr(PDESpecifier,"Steady")) {
			double const *const ADV_b = DB.ADV_b;
			for (size_t dim = 0; dim < d; dim++) {
			for (size_t n = 0; n < Nn; n++) {
				b[n+dim*Nn] = ADV_b[dim];
			}}
		}
	} else {
		EXIT_UNSUPPORTED;
		printf("%p\n",XYZ);
	}

	return b;
}
