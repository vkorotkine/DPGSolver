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

#ifndef DPG__numerical_flux_h__INCLUDED
#define DPG__numerical_flux_h__INCLUDED
/** \file
 *  \brief Provides containers and functions relating to the supported numerical fluxes.
 */

struct Numerical_Flux_Input;
struct mutable_Numerical_Flux;
struct Simulation;

#include <stdbool.h>
#include "boundary.h"

///\{ \name The maximum number of outputs from the numerical flux functions.
#define MAX_NUM_FLUX_OUT 3 ///< See the members of \ref Numerical_Flux.
///\}

/** \brief Function pointer to \ref Numerical_Flux computing functions.
 *  \return Standard.
 *
 *  \param num_flux_i \ref Numerical_Flux_Input.
 *  \param num_flux   \ref Numerical_Flux.
 */
typedef void (*compute_Numerical_Flux_fptr)
	(const struct Numerical_Flux_Input* num_flux_i,
	 struct mutable_Numerical_Flux* num_flux
	);

/// \brief Container holding data used for computing the numerical fluxes and their Jacobians.
struct Numerical_Flux_Input {
	struct Boundary_Value_Input bv_l; ///< \ref Boundary_Value_Input container.
	struct Boundary_Value       bv_r; ///< \ref Boundary_Value container.

	const bool has_1st_order, ///< \ref Test_Case::has_1st_order.
	           has_2nd_order; ///< \ref Test_Case::has_2nd_order.

	/// \ref compute_Numerical_Flux_fptr calling appropriate 1st/2nd order functions.
	compute_Numerical_Flux_fptr compute_Numerical_Flux;

	/// \ref compute_Numerical_Flux_fptr for the 1st order fluxes.
	compute_Numerical_Flux_fptr compute_Numerical_Flux_1st;

	/// \ref compute_Numerical_Flux_fptr for the 2nd order fluxes.
	compute_Numerical_Flux_fptr compute_Numerical_Flux_2nd;
};

/** \brief Container storing the numerical fluxes and their Jacobians.
 *  The number of members should be equal to \ref MAX_NUM_FLUX_OUT. */
struct Numerical_Flux {
	const struct const_Multiarray_d* nnf; ///< The normal numerical fluxes.

	/// \brief Container for information from either side of the face.
	struct Neigh_Info_NF {
		/// The Jacobian of the normal numerical fluxes wrt the solution variables.
		const struct const_Multiarray_d* dnnf_ds;

		/// The Jacobian of the normal numerical fluxes wrt the solution gradient variables.
		const struct const_Multiarray_d* dnnf_dg;
	} neigh_info[2]; ///< \ref Neigh_Info_NF.
};

/// \brief `mutable` version of \ref Numerical_Flux.
struct mutable_Numerical_Flux {
	struct Multiarray_d* nnf; ///< Defined in \ref Numerical_Flux.

	/// Defined in \ref Numerical_Flux.
	struct m_Neigh_Info_NF {
		struct Multiarray_d* dnnf_ds; ///< Defined in \ref Numerical_Flux.
		struct Multiarray_d* dnnf_dg; ///< Defined in \ref Numerical_Flux.
	} neigh_info[2]; ///< Defined in \ref Numerical_Flux.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Numerical_Flux_Input container.
 *  \return See brief. */
struct Numerical_Flux_Input* constructor_Numerical_Flux_Input
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Numerical_Flux_Input container.
void destructor_Numerical_Flux_Input
	(struct Numerical_Flux_Input* num_flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Numerical_Flux container.
 *  \return See brief. */
struct Numerical_Flux* constructor_Numerical_Flux
	(const struct Numerical_Flux_Input* num_flux_i ///< \ref Numerical_Flux_Input.
	);

/// \brief Destructor for a \ref Numerical_Flux container.
void destructor_Numerical_Flux
	(struct Numerical_Flux* num_flux ///< Standard.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr for 1st order numerical fluxes only.
void compute_Numerical_Flux_1
	(const struct Numerical_Flux_Input* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux* num_flux        ///< See brief.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr for 1st and 2nd order numerical fluxes.
void compute_Numerical_Flux_12
	(const struct Numerical_Flux_Input* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux* num_flux        ///< See brief.
	);

#endif // DPG__numerical_flux_h__INCLUDED
