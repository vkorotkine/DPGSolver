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

#include <assert.h>

#include "def_templates_vector.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void invert_Vector_T (struct Vector_T* a)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		assert(a->data[i] != 0.0);
		a->data[i] = 1.0/(a->data[i]);
	}
}

void add_to_Vector_T_T (struct Vector_T* a, const Type* b)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] += b[i];
}

void add_to_Vector_T (struct Vector_T* a, const struct const_Vector_T* b)
{
	assert(a->ext_0 == b->ext_0);
	const ptrdiff_t ext_0 = a->ext_0;
	for (int i = 0; i < ext_0; ++i)
		a->data[i] += b->data[i];
}

Type dot_product_Vector_T (const Type alpha, const struct const_Vector_T*const a, const struct const_Vector_T*const b)
{
	assert(a->ext_0 == b->ext_0);

	Type out = 0.0;
	for (int i = 0; i < a->ext_0; ++i)
		out += a->data[i]*b->data[i];

	return alpha*out;
}

Type dot_product_Vector_RT (const Type alpha, const struct const_Vector_R*const a, const struct const_Vector_T*const b)
{
	assert(a->ext_0 == b->ext_0);

	Type out = 0.0;
	for (int i = 0; i < a->ext_0; ++i)
		out += a->data[i]*b->data[i];

	return alpha*out;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
