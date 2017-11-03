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

// Interface functions ********************************************************************************************** //

void print_Matrix_c (const struct Matrix_c*const a)
{
	const double tol = EPS;

	const int n_dec = 4; // Number of places after the decimal.

	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	if (ext_0 == 0 || ext_1 == 0) {
		printf("Called print_Matrix_c for input with extents: [%td,%td].\n\n",ext_0,ext_1);
		return;
	}

	const double complex* data = a->data;

	char format_d[10];
	char format_i[10];
	sprintf(format_d,"%s%de%s","% .",n_dec," ");
	sprintf(format_i,"%s%dd%s","% ",n_dec+7," ");

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const double complex val = *data++;
				const double val_r = creal(val),
				             val_i = cimag(val);
//				printf("% .4e%c%.4ei ",Ar,(Ac >= 0.0)?'+':'\0',Ac);
				if (isnan(val_r) || (fabs(val_r) > tol)) {
					printf(format_d,val);
					printf("%c",(var_i > 0.0 ? '+':'\0'));
					printf(format_d,val);
					printf("i");
				}
				EXIT_UNSUPPORTED;
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
				const double complex val = data[i+ext_0*j];
				if (isnan(val) || (cabs(val) > tol))
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

void print_const_Matrix_c (const struct const_Matrix_c*const a)
{
	print_Matrix_c((struct Matrix_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
