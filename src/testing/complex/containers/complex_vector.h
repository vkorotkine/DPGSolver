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

#ifndef DPG__complex_vector_h__INCLUDED
#define DPG__complex_vector_h__INCLUDED
/** \file
 *  \brief Provides `complex` Vector_\* containers and related functions.
 *
 *  Potentially relevant comments may be found in \ref vector.h.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_vector_constructors.h"
//#include "complex_vector_math.h"
//#include "complex_vector_print.h"

/// \brief `complex` version of \ref Vector_d.
struct Vector_c {
	ptrdiff_t ext_0; ///< See brief.

	bool owns_data;       ///< See brief.
	double complex* data; ///< See brief.
};

/// \brief `const` version of \ref Vector_c.
struct const_Vector_c {
	const ptrdiff_t ext_0; ///< See brief.

	const bool owns_data;            ///< See brief.
	const double complex*const data; ///< See brief.
};

// Interface functions ********************************************************************************************** //

#endif // DPG__complex_vector_h__INCLUDED
