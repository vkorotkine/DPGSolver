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

/** \brief `complex` version of \ref get_col_Matrix_d.
 *  \return See brief. */
double complex* get_col_Matrix_c
	(const ptrdiff_t col,     ///< See brief.
	 const struct Matrix_c* a ///< See brief.
	);

#endif // DPG__complex_matrix_h__INCLUDED
