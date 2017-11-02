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

#ifndef DPG__complex_multiarray_h__INCLUDED
#define DPG__complex_multiarray_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* containers and related functions for complex variables.
 *
 *  See \ref multiarray.h for potentially relevant comments.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

#include "complex_multiarray_constructors.h"

struct const_Multiarray_d;

/// \brief Multiarray (`double complex`).
struct Multiarray_c {
	char layout; ///< Defined in \ref Multiarray_d.

	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data; ///< Defined in \ref Multiarray_d.

	double complex* data; ///< Defined in \ref Multiarray_d.
};

/// \brief `const` version of \ref Multiarray_c.
struct const_Multiarray_c {
	const char layout; ///< Defined in \ref Multiarray_c.

	const int order;               ///< Defined in \ref Multiarray_c.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_c.

	const bool owns_data;            ///< Defined in \ref Multiarray_c.
	const double complex*const data; ///< Defined in \ref Multiarray_c.
};

// Interface functions ********************************************************************************************** //

/// \brief Set the data of the \ref Multiarary_c\* container from that of the \ref Multiarray_d\* container.
void set_Multiarray_c_Multiarray_d
	(struct Multiarray_c* a,            ///< Multiarray with data to be set.
	 const struct const_Multiarray_d* b ///< Multiarray from which to take data.
	);

#endif // DPG__complex_multiarray_h__INCLUDED
