// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix_print.h"

#include <stdio.h>
#include <math.h>

#include "macros.h"

#include "matrix.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Matrix_d (const struct Matrix_d*const a, const double tol)
{
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	const double* data = a->data;

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double val = *data++;
				printf("% .12e ",( (isnan(val) || (fabs(val) > tol)) ? val : 0.0 ));
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double val = data[i+ext_0*j];
				printf("% .4e ",( (isnan(val) || (fabs(val) > tol)) ? val : 0.0 ));
			}
			printf("\n");
		}
		printf("\n");
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void print_const_Matrix_d (const struct const_Matrix_d*const a, const double tol)
{
	struct Matrix_d* local = constructor_move_Matrix_d_d(a->layout,a->ext_0,a->ext_1,false,(double*)a->data);
	print_Matrix_d(local,tol);
	free(local);
}

void print_Matrix_i (const struct Matrix_i*const a)
{
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	const int* data = a->data;

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const int val = *data++;
				printf("% 12d ",val);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const int val = data[i+ext_0*j];
				printf("% 12d ",val);
			}
			printf("\n");
		}
		printf("\n");
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void print_const_Matrix_i (const struct const_Matrix_i*const a)
{
	struct Matrix_i* local = constructor_move_Matrix_i_i(a->layout,a->ext_0,a->ext_1,false,(int*)a->data);
	print_Matrix_i(local);
	free(local);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
