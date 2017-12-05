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

#ifndef DPG__complex_multiarray_minimal_h__INCLUDED
#define DPG__complex_multiarray_minimal_h__INCLUDED
/** \file
 *  \brief Provides **minimal** \ref Multiarray_c containers for complex variables.
 *
 *  See \ref multiarray.h for potentially relevant comments.
 *
 *  As these containers are only used for testing purposes, the minimal amount of code is exposed here as required by
 *  functions in the core of the src (i.e. all directories exclusing src/testing). All other required functions for
 *  testing are implemented in the corresponding test_support_\* file(s).
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_multiarray_minimal_constructors.h"

struct const_Multiarray_d;

/// \brief Multiarray (`double complex`).
struct Multiarray_c {
	char layout; ///< Defined in \ref Multiarray_T.

	int order;          ///< Defined in \ref Multiarray_T.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_T.

	bool owns_data; ///< Defined in \ref Multiarray_T.

	double complex* data; ///< Defined in \ref Multiarray_T.
};

/// \brief `const` version of \ref Multiarray_c.
struct const_Multiarray_c {
	const char layout; ///< Defined in \ref Multiarray_c.

	const int order;               ///< Defined in \ref Multiarray_c.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_c.

	const bool owns_data;            ///< Defined in \ref Multiarray_c.
	const double complex*const data; ///< Defined in \ref Multiarray_c.
};

#endif // DPG__complex_multiarray_minimal_h__INCLUDED
