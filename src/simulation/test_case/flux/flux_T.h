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
 *  \brief Provides templated containers and functions relating to fluxes of the supported PDEs.
 *
 *  The memory layout of the fluxes (node, dimension, equation) was chosen such that the memory stride is minimized when
 *  converting from physical to reference space.
 */

#include "def_templates_multiarray.h"
#include "def_templates_flux.h"

struct Flux_Input_T;
struct mutable_Flux_T;
struct Simulation;

#include <stdbool.h>

/** \brief Function pointer to \ref Flux_T computing functions.
 *  \return Standard.
 *
 *  \param flux_i \ref Flux_Input_T.
 *  \param flux   \ref Flux_T.
 */
typedef void (*compute_Flux_fptr_T)
	(const struct Flux_Input_T* flux_i,
	 struct mutable_Flux_T* flux
	);

/// \brief Container holding data used for computing the fluxes and flux Jacobians.
struct Flux_Input_T {
	/** Array of flags for which of the \ref Flux_T members should be computed. The entries correspond to the
	 *  following members in \ref Flux_T.
	 *  - [0]: f
	 *  - [1]: df_ds
	 *  - [2]: df_dg
	 *  - [3]: d2f_ds2
	 *  - [4]: d2f_dsdg
	 *  - [5]: d2f_dg2
	 */
	const bool* compute_member;

	const int n_eq,  ///< \ref Test_Case_T::n_eq.
	          n_var; ///< \ref Test_Case_T::n_var.

	const bool has_1st_order, ///< \ref Test_Case_T::has_1st_order.
	           has_2nd_order; ///< \ref Test_Case_T::has_2nd_order.

	const struct const_Multiarray_T* s;   ///< The solution variables.
	const struct const_Multiarray_T* g;   ///< The solution gradient variables.
	const struct const_Multiarray_T* xyz; ///< The xyz-coordinates.

	compute_Flux_fptr_T compute_Flux;     ///< \ref compute_Flux_fptr_T calling appropriate 1st/2nd order functions.
	compute_Flux_fptr_T compute_Flux_1st; ///< \ref compute_Flux_fptr_T for the 1st order fluxes.
	compute_Flux_fptr_T compute_Flux_2nd; ///< \ref compute_Flux_fptr_T for the 2nd order fluxes.
};

/** \brief Container storing the fluxes and flux Jacobians.
 *  The number of members should be equal to \ref MAX_FLUX_OUT. */
struct Flux_T {
	const struct const_Multiarray_T* f;       ///< The fluxes.
	const struct const_Multiarray_T* df_ds;   ///< The Jacobian of the fluxes wrt the solution variables.
	const struct const_Multiarray_T* df_dg;   ///< The Jacobian of the fluxes wrt the solution gradient variables.
	const struct const_Multiarray_T* d2f_ds2; ///< The Hessian of the fluxes wrt the solution variables.
};

/// \brief `mutable` version of \ref Flux_T.
struct mutable_Flux_T {
	struct Multiarray_T* f;       ///< See brief.
	struct Multiarray_T* df_ds;   ///< See brief.
	struct Multiarray_T* df_dg;   ///< See brief.
	struct Multiarray_T* d2f_ds2; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Flux_Input_T container.
 *  \return See brief. */
struct Flux_Input_T* constructor_Flux_Input_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Flux_Input_T using \ref Test_Case_T::solver_method_curr = 'e'.
 *  \return See brief. */
struct Flux_Input_T* constructor_Flux_Input_T_e
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Flux_Input_T container.
void destructor_Flux_Input_T
	(struct Flux_Input_T* flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Flux_T container.
 *  \return See brief. */
struct Flux_T* constructor_Flux_T
	(const struct Flux_Input_T* flux_i ///< \ref Flux_Input_T.
	);

/// \brief Destructor for a \ref Flux_T container.
void destructor_Flux_T
	(struct Flux_T* flux ///< Standard.
		);

/// \brief Destructor for a \ref Flux_T container if the input is non-NULL, otherwise simply return.
void destructor_conditional_Flux_T
	(struct Flux_T* flux ///< Standard.
	);

/// \brief Version of \ref compute_Flux_fptr_T for 1st order fluxes only.
void compute_Flux_1_T
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);

/// \brief Version of \ref compute_Flux_fptr_T for 2nd order fluxes only.
void compute_Flux_2_T
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);

/// \brief Version of \ref compute_Flux_fptr_T for 1st and 2nd order fluxes.
void compute_Flux_12_T
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);


#include "undef_templates_multiarray.h"
#include "undef_templates_flux.h"
