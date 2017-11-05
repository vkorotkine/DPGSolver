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
	char layout; ///< Defined in \ref Matrix_d.

	ptrdiff_t ext_0, ///< Defined in \ref Matrix_d.
	          ext_1; ///< Defined in \ref Matrix_d.

	bool owns_data; ///< Defined in \ref Matrix_d.
	double complex* data; ///< Defined in \ref Matrix_d.
};

/// \brief `const` version of \ref Matrix_c.
struct const_Matrix_c {
	const char layout; ///< Defined in \ref Matrix_c.

	const ptrdiff_t ext_0, ///< Defined in \ref Matrix_c.
	                ext_1; ///< Defined in \ref Matrix_c.

	const bool owns_data;    ///< Defined in \ref Matrix_c.
	const double complex*const data; ///< Defined in \ref Matrix_c.
};

// Interface functions ********************************************************************************************** //

#endif // DPG__complex_matrix_h__INCLUDED
