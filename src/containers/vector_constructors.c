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

#include "vector_constructors.h"

#include "macros.h"
#include "definitions_tol.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "file_processing.h"
#include "math_functions.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "vector_constructors_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "def_templates_type_i.h"
#include "def_templates_matrix_i.h"
#include "def_templates_multiarray_i.h"
#include "def_templates_vector_i.h"
#include "vector_constructors_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

const struct const_Vector_i* constructor_const_Vector_i_inds_eq_1_const_Matrix_d (const struct const_Matrix_d*const a)
{
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	struct Vector_i*const inds = constructor_empty_Vector_i(ext_0); // returned
	for (int i = 0; i < ext_0; ++i) {
		bool found = 0;
		const double*const op_row = get_row_const_Matrix_d(i,a);
		for (int j = 0; j < ext_1; ++j) {
			if (equal_d(op_row[j],1.0,4e0*EPS)) {
				found = 1;
				inds->data[i] = j;
				break;
			}
		}
		if (!found) {
			printf("Did not find the column corresponding to the 1.0 for row: %d.\n",i);
			print_const_Matrix_d(a);
			for (int j = 0; j < ext_1; ++j)
				printf("%d % .3e % .3e\n",j,op_row[j],op_row[j]-1.0);
			EXIT_UNSUPPORTED;
		}
	}
	return (struct const_Vector_i*) inds;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
