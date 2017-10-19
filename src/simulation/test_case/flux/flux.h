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

#ifndef DPG__flux_h__INCLUDED
#define DPG__flux_h__INCLUDED
/** \file
 *  \brief Provides containers and functions relating to fluxes of the supported PDEs.
 */

struct Flux_Input;
struct mutable_Flux;
struct Simulation;

#include <stdbool.h>
#include "flux_euler.h"

///\{ \name The maximum number of outputs from the flux functions.
#define MAX_FLUX_OUT 3 ///< See the members of \ref Flux.
///\}

/** \brief Function pointer to \ref Flux computing functions.
 *  \return Standard.
 *
 *  \param flux_i \ref Flux_Input.
 *  \param flux   \ref Flux.
 */
typedef void (*compute_Flux_fptr)
	(const struct Flux_Input* flux_i,
	 struct mutable_Flux* flux
	);

/// \brief Container holding data used for computing the fluxes and flux Jacobians.
struct Flux_Input {
	const bool* compute_member; ///< Array of flags for which of the \ref Flux members should be computed.

	const int d,     ///< \ref Simulation::d.
	          n_eq,  ///< \ref Test_Case::n_eq.
	          n_var; ///< \ref Test_Case::n_var.

	const bool has_1st_order, ///< \ref Test_Case::has_1st_order.
	           has_2nd_order; ///< \ref Test_Case::has_2nd_order.

	const struct const_Multiarray_d* s;   ///< The solution variables.
	const struct const_Multiarray_d* g;   ///< The solution gradient variables.
	const struct const_Multiarray_d* xyz; ///< The xyz-coordinates.

	compute_Flux_fptr compute_Flux;     ///< \ref compute_Flux_fptr calling appropriate 1st/2nd order functions.
	compute_Flux_fptr compute_Flux_1st; ///< \ref compute_Flux_fptr for the 1st order fluxes.
	compute_Flux_fptr compute_Flux_2nd; ///< \ref compute_Flux_fptr for the 2nd order fluxes.
};

/** \brief Container storing the fluxes and flux Jacobians.
 *  The number of members should be equal to \ref MAX_FLUX_OUT. */
struct Flux {
	const struct const_Multiarray_d* f;     ///< The fluxes.
	const struct const_Multiarray_d* df_ds; ///< The Jacobian of the fluxes wrt the solution variables.
	const struct const_Multiarray_d* df_dg; ///< The Jacobian of the fluxes wrt the solution gradient variables.
};

/// \brief `mutable` version of \ref Flux.
struct mutable_Flux {
	struct Multiarray_d* f;     ///< Defined in \ref Flux.
	struct Multiarray_d* df_ds; ///< Defined in \ref Flux.
	struct Multiarray_d* df_dg; ///< Defined in \ref Flux.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Flux_Input container.
 *  \return See brief. */
struct Flux_Input* constructor_Flux_Input
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Flux_Input container.
void destructor_Flux_Input
	(struct Flux_Input* flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Flux container.
 *  \return See brief. */
struct Flux* constructor_Flux
	(const struct Flux_Input* flux_i ///< \ref Flux_Input.
	);

/// \brief Destructor for a \ref Flux container.
void destructor_Flux
	(struct Flux* flux ///< Standard.
	);

/// \brief Version of \ref compute_Flux_fptr for 1st order fluxes only.
void compute_Flux_1
	(const struct Flux_Input* flux_i, ///< Defined for \ref compute_Flux_fptr.
	 struct mutable_Flux* flux        ///< Defined for \ref compute_Flux_fptr.
	);

/// \brief Version of \ref compute_Flux_fptr for 1st and 2nd order fluxes.
void compute_Flux_12
	(const struct Flux_Input* flux_i, ///< Defined for \ref compute_Flux_fptr.
	 struct mutable_Flux* flux        ///< Defined for \ref compute_Flux_fptr.
	);

#endif // DPG__flux_h__INCLUDED
