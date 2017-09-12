// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "vector_print.h"

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

void print_Vector_d (const struct Vector_d*const a, const double tol)
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

void print_const_Vector_d (const struct const_Vector_d*const a, const double tol)
{
	print_Vector_d((struct Vector_d*)a,tol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
