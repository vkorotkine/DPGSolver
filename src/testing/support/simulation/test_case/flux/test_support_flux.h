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

#ifndef DPG__test_support_flux_h__INCLUDED
#define DPG__test_support_flux_h__INCLUDED
/** \file
 *  \brief Provides `complex` versions of containers and functions defined in \ref flux.h.
 *
 *  As the complex functions are never used to compute the linearization terms, they are omitted from the containers.
 */

struct Flux_Input_c;
struct mutable_Flux_c;
struct Simulation;

#include "flux.h"

/** \brief `complex` version of \ref compute_Flux_fptr.
 *  \return Standard.
 *
 *  \param flux_i Defined for \ref compute_Flux_fptr.
 *  \param flux   Defined for \ref compute_Flux_fptr.
 */
typedef void (*compute_Flux_c_fptr)
	(const struct Flux_Input_c* flux_i,
	 struct mutable_Flux_c* flux
	);

/// \brief Derived `complex` version of \ref Flux_Input.
struct Flux_Input_c {
	struct Flux_Input flux_i; ///< Base \ref Flux_Input.

	const struct const_Multiarray_c* s;   ///< Defined for \ref Flux_Input.
	const struct const_Multiarray_c* g;   ///< Defined for \ref Flux_Input.
	const struct const_Multiarray_c* xyz; ///< Defined for \ref Flux_Input.

	compute_Flux_c_fptr compute_Flux;     ///< Defined for \ref Flux_Input.
	compute_Flux_c_fptr compute_Flux_1st; ///< Defined for \ref Flux_Input.
	compute_Flux_c_fptr compute_Flux_2nd; ///< Defined for \ref Flux_Input.
};

/// \brief `complex` version of \ref Flux.
struct Flux_c {
	const struct const_Multiarray_c* f; ///< Defined for \ref Flux.
};

/// \brief `mutable` version of \ref Flux_c.
struct mutable_Flux_c {
	struct Multiarray_c* f; ///< Defined in \ref Flux_c.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Flux_Input_c container.
 *  \return See brief. */
struct Flux_Input_c* constructor_Flux_Input_c
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Flux_Input_c container.
void destructor_Flux_Input_c
	(struct Flux_Input_c* flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Flux_c container.
 *  \return See brief. */
struct Flux_c* constructor_Flux_c
	(const struct Flux_Input_c* flux_i ///< \ref Flux_Input.
	);

/// \brief Destructor for a \ref Flux_c container.
void destructor_Flux_c
	(struct Flux_c* flux ///< Standard.
	);

/// \brief `complex` version of \ref compute_Flux_1.
void compute_Flux_c_1
	(const struct Flux_Input_c* flux_i, ///< Defined for \ref compute_Flux_c_fptr.
	 struct mutable_Flux_c* flux        ///< Defined for \ref compute_Flux_c_fptr.
	);

/// \brief `complex` version of \ref compute_Flux_12.
void compute_Flux_c_12
	(const struct Flux_Input_c* flux_i, ///< Defined for \ref compute_Flux_c_fptr.
	 struct mutable_Flux_c* flux        ///< Defined for \ref compute_Flux_c_fptr.
	);

#endif // DPG__test_support_flux_h__INCLUDED
