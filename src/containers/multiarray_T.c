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

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

/** \brief Get pointer to col of col-major Multiarray_\*.
 *  \return See brief. */
Type* get_col_Multiarray_T
	(const ptrdiff_t col,   ///< Desired column.
	 struct Multiarray_T* a ///< Multiarray.
	)
{
	assert(a->layout == 'C');

	const ptrdiff_t ext_0 = a->extents[0];
	return &a->data[col*ext_0];
}

/** \brief `const` version of \ref get_col_Multiarray_T.
 *  \return See brief. */
const Type* get_col_const_Multiarray_T
	(const ptrdiff_t col,               ///< See brief.
	 const struct const_Multiarray_T* a ///< See brief.
	)
{
	return (const Type*) get_col_Multiarray_T(col,(struct Multiarray_T*)a);
}

/// \brief Set all data entries to the input value.
void set_to_value_Multiarray_T
	(struct Multiarray_T*const a, ///< Multiarray.
	 const Type val               ///< Value.
	)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		a->data[i] = val;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
