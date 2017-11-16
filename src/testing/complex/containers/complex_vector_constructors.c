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

#include "complex_vector_constructors.h"

#include <assert.h>
#include <stdlib.h>

#include "complex_vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Default constructors ********************************************************************************************* //

// Empty constructors *********************************************************************************************** //

// Zero constructors ************************************************************************************************ //

struct Vector_c* constructor_zero_Vector_c (const ptrdiff_t ext_0)
{
	double complex* data = calloc(ext_0 , sizeof *data); // keep
	return constructor_move_Vector_c_c(ext_0,true,data);
}

// Copy constructors ************************************************************************************************ //

// Move constructors ************************************************************************************************ //

struct Vector_c* constructor_move_Vector_c_c (const ptrdiff_t ext_0, const bool owns_data, double complex*const data)
{
	struct Vector_c* dest = calloc(1,sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

// Set constructors ************************************************************************************************* //

// Special constructors ********************************************************************************************* //

// Destructors ****************************************************************************************************** //

void destructor_Vector_c (struct Vector_c* a)
{
	assert(a != NULL);

	if (a->owns_data)
		free(a->data);
	free(a);
}

void destructor_const_Vector_c (const struct const_Vector_c* a)
{
	destructor_Vector_c((struct Vector_c*)a);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
