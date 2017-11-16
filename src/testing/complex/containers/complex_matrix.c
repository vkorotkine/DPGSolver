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

#include "complex_matrix.h"

#include <assert.h>

#include "macros.h"

// Static function declarations ************************************************************************************* //

/** \brief `complex` version of \ref set_value_fptr.
 *
 *  \param dest See brief.
 *  \param src  See brief.
 */
typedef void (*set_value_c_fptr)
	(double complex*const dest,
	 const double complex src
	);

/// \brief `complex` version of \ref set_value_insert.
void set_value_insert_c
	(double complex*const dest, ///< See brief.
	 const double complex src   ///< See brief.
	);

/// \brief `complex` version of \ref set_value_add.
void set_value_add_c
	(double complex*const dest, ///< See brief.
	 const double complex src   ///< See brief.
	);

// Interface functions ********************************************************************************************** //

double complex* get_row_Matrix_c (const ptrdiff_t row, const struct Matrix_c* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

const double complex* get_row_const_Matrix_c (const ptrdiff_t row, const struct const_Matrix_c* a)
{
	assert(a->layout == 'R');
	return &a->data[row*(a->ext_1)];
}

double complex* get_col_Matrix_c (const ptrdiff_t col, const struct Matrix_c* a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

const double complex* get_col_const_Matrix_c (const ptrdiff_t col, const struct const_Matrix_c* a)
{
	assert(a->layout == 'C');
	return &a->data[col*(a->ext_0)];
}

void set_to_value_Matrix_c (struct Matrix_c*const a, const double complex val)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

void set_block_Matrix_c
	(struct Matrix_c* a, const struct const_Matrix_c* a_sub, const ptrdiff_t row0, const ptrdiff_t col0,
	 const char set_type)
{
	assert(a->layout == a_sub->layout); // Add support if required.
	assert(a->layout == 'R');

	set_value_c_fptr set_value = NULL;
	switch (set_type) {
		case 'i': set_value = set_value_insert_c; break;
		case 'a': set_value = set_value_add_c;    break;
		default:  EXIT_ERROR("Unsupported: %c.\n",set_type); break;
	}

	const ptrdiff_t ext_0 = a_sub->ext_0,
	                ext_1 = a_sub->ext_1;

	assert(row0+ext_0 <= a->ext_0);
	assert(col0+ext_1 <= a->ext_1);

	for (int i = 0, row = row0; i < ext_0; ++i, ++row) {
		const double complex*const data_as = get_row_const_Matrix_c(i,a_sub);

		double complex* data_a = get_row_Matrix_c(row,a);
		data_a += col0;
		for (int j = 0; j < ext_1; ++j)
			set_value(&data_a[j],data_as[j]);
	}
}

void set_block_Matrix_c_d
	(struct Matrix_c* a, const struct const_Matrix_d* a_sub, const ptrdiff_t row0, const ptrdiff_t col0,
	 const char set_type)
{
	const struct const_Matrix_c*const a_sub_c = constructor_copy_const_Matrix_c_Matrix_d(a_sub); // destructed
	set_block_Matrix_c(a,a_sub_c,row0,col0,set_type);
	destructor_const_Matrix_c(a_sub_c);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

void set_value_insert_c (double complex*const dest, const double complex src)
{
	*dest = src;
}

void set_value_add_c (double complex*const dest, const double complex src)
{
	*dest += src;
}
