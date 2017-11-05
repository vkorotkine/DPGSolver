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

#include "complex_matrix_print.h"

#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"

#include "complex_matrix.h"

// Static function declarations ************************************************************************************* //

///\{ \name The number of decimal places to display.
#define N_DEC 4
///\}

/// \brief Print a real value to the terminal with default format.
static void print_real
	(const double complex val ///< The complex value.
	);

/// \brief Print an imaginary value to the terminal with default format.
static void print_imag
	(const double complex val ///< The complex value.
	);

// Interface functions ********************************************************************************************** //

void print_Matrix_c (const struct Matrix_c*const a)
{

	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	if (ext_0 == 0 || ext_1 == 0) {
		printf("Called print_Matrix_c for input with extents: [%td,%td].\n\n",ext_0,ext_1);
		return;
	}

	const double complex* data = a->data;
	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double complex val = *data++;
				print_real(val);
				print_imag(val);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double complex val = data[i+ext_0*j];
				print_real(val);
				print_imag(val);
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

void print_const_Matrix_c (const struct const_Matrix_c*const a)
{
	print_Matrix_c((struct Matrix_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void print_real (const double complex val)
{
	static char format_d[10],
	            format_i[10];

	static bool format_set = false;
	if (!format_set) {
		format_set = true;
		sprintf(format_d,"%s%de","% .",N_DEC);
		sprintf(format_i,"%s%dd","% ",N_DEC+7);
	}

	static const double tol = EPS;
	const double val_r = creal(val);
	if (isnan(val_r) || (fabs(val_r) > tol))
		printf(format_d,val_r);
	else
		printf(format_i,0);
}

static void print_imag (const double complex val)
{
	static char format_dp[10],
	            format_dm[10],
	            format_ip[10];

	static bool format_set = false;
	if (!format_set) {
		format_set = true;
		sprintf(format_dp,"%s%de%s","+%.",N_DEC,"i ");
		sprintf(format_dm,"%s%de%s","% .",N_DEC,"i ");
		sprintf(format_ip,"%s%dd%s","+% ",N_DEC+6,"i ");
	}

	static const double tol_i = EPS*EPS*EPS;
	const double val_i = cimag(val);
	if (isnan(val_i) || (fabs(val_i) > tol_i)) {
		if (val_i >= 0.0)
			printf(format_dp,val_i);
		else
			printf(format_dm,val_i);
	} else {
		printf(format_ip,0);
	}
}
