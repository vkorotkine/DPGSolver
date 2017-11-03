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

#include "complex_matrix_constructors.h"

#include <stdlib.h>
#include <assert.h>

#include "complex_matrix.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

/** \brief Copy constructor for a \ref Matrix_c\* from a \ref Matrix_d\*.
 *  \return See brief. */
static struct Matrix_c* constructor_copy_Matrix_c_Matrix_d
	(struct Matrix_d* src ///< The source matrix.
	);

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

// Copy constructors ************************************************************************************************ //

const struct const_Matrix_c* constructor_copy_const_Matrix_c_Matrix_d (const struct const_Matrix_d* src)
{
	return (struct const_Matrix_c*) constructor_copy_Matrix_c_Matrix_d((struct Matrix_d*)src);
}

// Move constructors ************************************************************************************************ //

struct Matrix_c* constructor_move_Matrix_c_c
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, double complex*const data)
{
    struct Matrix_c* dest = calloc(1,sizeof *dest); // returned

    dest->layout    = layout;
    dest->ext_0     = ext_0;
    dest->ext_1     = ext_1;
    dest->owns_data = owns_data;
    dest->data      = data;

    return dest;
}

const struct const_Matrix_c* constructor_move_const_Matrix_c_c
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data,
	 const double complex*const data)
{
	return (const struct const_Matrix_c*)
		constructor_move_Matrix_c_c(layout,ext_0,ext_1,owns_data,(double complex*)data);
}

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

void destructor_Matrix_c (struct Matrix_c* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Matrix_c (const struct const_Matrix_c* a)
{
	destructor_Matrix_c((struct Matrix_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Matrix_c* constructor_copy_Matrix_c_Matrix_d (struct Matrix_d* src)
{
	const ptrdiff_t size = (src->ext_0)*(src->ext_1);
	const double*const data_src = src->data;

	double complex* data = calloc(size , sizeof *data); // keep
//	double complex* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < size; i++)
		data[i] = data_src[i];

	return constructor_move_Matrix_c_c(src->layout,src->ext_0,src->ext_1,true,data);
}
