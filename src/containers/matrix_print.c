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
/** \file
 */

#include "matrix_print.h"

#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"

#include "matrix.h"

// Static function declarations ************************************************************************************* //

/** \brief Check whether either of the matrix extents is zero.
 *  \return `true` if yes; `false` otherwise. */
static bool check_Matrix_extents_zero
	(const ptrdiff_t ext_0, ///< \ref First extent.
	 const ptrdiff_t ext_1  ///< \ref Second extent.
	);

// Interface functions ********************************************************************************************** //

void print_Matrix_d_tol (const struct Matrix_d*const a, const double tol)
{
	const int n_dec = 4; // Number of places after the decimal.

	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	if (check_Matrix_extents_zero(ext_0,ext_1))
		printf("Called print_Matrix_d for input with extents: [%td,%td].\n\n",ext_0,ext_1);

	const double* data = a->data;

	char format_d[10];
	char format_i[10];
	sprintf(format_d,"%s%de%s","% .",n_dec," ");
	sprintf(format_i,"%s%dd%s","% ",n_dec+7," ");

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double val = *data++;
				if (isnan(val) || (fabs(val) > tol))
					printf(format_d,val);
				else
					printf(format_i,0);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double val = data[i+ext_0*j];
				if (isnan(val) || (fabs(val) > tol))
					printf(format_d,val);
				else
					printf(format_i,0);
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

void print_const_Matrix_d_tol (const struct const_Matrix_d*const a, const double tol)
{
	print_Matrix_d_tol((struct Matrix_d*)a,tol);
}

void print_Matrix_d (const struct Matrix_d*const a)
{
	print_Matrix_d_tol(a,EPS);
}

void print_const_Matrix_d (const struct const_Matrix_d*const a)
{
	print_Matrix_d((const struct Matrix_d*)a);
}

void print_Matrix_i (const struct Matrix_i*const a)
{
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	if (check_Matrix_extents_zero(ext_0,ext_1))
		printf("Called print_Matrix_d for input with extents: [%td,%td].\n\n",ext_0,ext_1);

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

static bool check_Matrix_extents_zero (const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	if (ext_0 == 0 || ext_1 == 0)
		return true;
	return false;
}
