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

#ifndef DPG__complex_matrix_h__INCLUDED
#define DPG__complex_matrix_h__INCLUDED
/** \file
 *  \brief Provides `complex` Matrix_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref matrix.h.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_matrix_constructors.h"
#include "complex_matrix_math.h"
#include "complex_matrix_print.h"

/// \brief `complex` version of \ref Matrix_d.
struct Matrix_c {
	char layout; ///< See brief.

	ptrdiff_t ext_0, ///< See brief.
	          ext_1; ///< See brief.

	bool owns_data;       ///< See brief.
	double complex* data; ///< See brief.
};

/// \brief `const` version of \ref Matrix_c.
struct const_Matrix_c {
	const char layout; ///< See brief.

	const ptrdiff_t ext_0, ///< See brief.
	                ext_1; ///< See brief.

	const bool owns_data;            ///< See brief.
	const double complex*const data; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief `complex` version of \ref get_row_Matrix_d.
 *  \return See brief. */
double complex* get_row_Matrix_c
	(const ptrdiff_t row,     ///< See brief.
	 const struct Matrix_c* a ///< See brief.
	);

/** \brief `complex` version of \ref get_row_const_Matrix_d.
 *  \return See brief. */
const double complex* get_row_const_Matrix_c
	(const ptrdiff_t row,           ///< See brief.
	 const struct const_Matrix_c* a ///< See brief.
	);

/** \brief `complex` version of \ref get_col_Matrix_d.
 *  \return See brief. */
double complex* get_col_Matrix_c
	(const ptrdiff_t col,     ///< See brief.
	 const struct Matrix_c* a ///< See brief.
	);

/** \brief `complex` version of \ref get_col_const_Matrix_d.
 *  \return See brief. */
const double complex* get_col_const_Matrix_c
	(const ptrdiff_t col,           ///< See brief.
	 const struct const_Matrix_c* a ///< See brief.
	);

/// \brief Set all data entries to the input value.
void set_to_value_Matrix_c
	(struct Matrix_c*const a, ///< Standard.
	 const double complex val ///< The value.
	);

/// \brief `complex` version of \ref set_block_Matrix_d.
void set_block_Matrix_c
	(struct Matrix_c* a,                 ///< See brief.
	 const struct const_Matrix_c* a_sub, ///< See brief.
	 const ptrdiff_t row0,               ///< See brief.
	 const ptrdiff_t col0,               ///< See brief.
	 const char set_type                 ///< See brief.
	);

/// \brief `complex` version of \ref set_block_Matrix_d with a `double` input block.
void set_block_Matrix_c_d
	(struct Matrix_c* a,                 ///< See brief.
	 const struct const_Matrix_d* a_sub, ///< See brief.
	 const ptrdiff_t row0,               ///< See brief.
	 const ptrdiff_t col0,               ///< See brief.
	 const char set_type                 ///< See brief.
	);

#endif // DPG__complex_matrix_h__INCLUDED
