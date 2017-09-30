/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/// \file

#include "vector_print.h"

#include "definitions_tol.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Vector_i (const struct Vector_i*const a)
{
	const ptrdiff_t ext = a->ext_0;

	const int* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
		printf("% 12d ",*data++);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_const_Vector_i (const struct const_Vector_i*const a)
{
	struct Vector_i* local = constructor_move_Vector_i_i(a->ext_0,false,(int*)a->data); // free
	print_Vector_i(local);
	free(local);
}

void print_Vector_d_tol (const struct Vector_d*const a, const double tol)
{
	const ptrdiff_t ext = a->ext_0;

	const double* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
		const double val = *data++;
		printf("% .4e ",( (isnan(val) || (fabs(val) > tol)) ? val : 0.0 ));
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_const_Vector_d_tol (const struct const_Vector_d*const a, const double tol)
{
	print_Vector_d_tol((struct Vector_d*)a,tol);
}

void print_Vector_d (const struct Vector_d*const a)
{
	print_Vector_d_tol(a,EPS);
}

void print_const_Vector_d (const struct const_Vector_d*const a)
{
	print_Vector_d((const struct Vector_d*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
