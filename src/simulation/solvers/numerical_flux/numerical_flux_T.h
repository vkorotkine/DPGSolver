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
 *  \brief Provides templated containers and functions relating to the supported numerical fluxes.
 */

struct Numerical_Flux_Input_T;
struct mutable_Numerical_Flux_T;
struct Simulation;

#include <stdbool.h>

/** \brief Function pointer to \ref Numerical_Flux_T computing functions.
 *  \return Standard.
 *
 *  \param num_flux_i \ref Numerical_Flux_Input_T.
 *  \param num_flux   \ref Numerical_Flux_T.
 */
typedef void (*compute_Numerical_Flux_fptr_T)
	(const struct Numerical_Flux_Input_T* num_flux_i,
	 struct mutable_Numerical_Flux_T* num_flux
	);

/// \brief Container holding data used for computing the numerical fluxes and their Jacobians.
struct Numerical_Flux_Input_T {
	struct Boundary_Value_Input_T bv_l; ///< \ref Boundary_Value_Input_T container.
	struct Boundary_Value_T       bv_r; ///< \ref Boundary_Value_T container.

	const int method; ///< \ref Simulation::method.

	const bool has_1st_order, ///< \ref Test_Case_T::has_1st_order.
	           has_2nd_order; ///< \ref Test_Case_T::has_2nd_order.

	/// \ref compute_Numerical_Flux_fptr_T calling appropriate 1st/2nd order functions.
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux;

	/// \ref compute_Numerical_Flux_fptr_T for the 1st order fluxes.
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux_1st;

	/// \ref compute_Numerical_Flux_fptr_T for the 2nd order fluxes.
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux_2nd;
};

/** \brief Container storing the numerical fluxes and their Jacobians.
 *  The number of members should be equal to \ref MAX_FLUX_OUT. */
struct Numerical_Flux_T {
	const struct const_Multiarray_T* nnf; ///< The normal numerical fluxes.

	/// \brief Container for information from either side of the face.
	struct Neigh_Info_NF_T {
		/// The Jacobian of the normal numerical fluxes wrt the solution variables.
		const struct const_Multiarray_T* dnnf_ds;

		/// The Jacobian of the normal numerical fluxes wrt the solution gradient variables.
		const struct const_Multiarray_T* dnnf_dg;
	} neigh_info[2]; ///< \ref Neigh_Info_NF_T.
};

/// \brief `mutable` version of \ref Numerical_Flux_T.
struct mutable_Numerical_Flux_T {
	struct Multiarray_T* nnf; ///< Defined in \ref Numerical_Flux_T.

	/// Defined in \ref Numerical_Flux_T.
	struct m_Neigh_Info_NF_T {
		struct Multiarray_T* dnnf_ds; ///< Defined in \ref Numerical_Flux_T.
		struct Multiarray_T* dnnf_dg; ///< Defined in \ref Numerical_Flux_T.
	} neigh_info[2]; ///< Defined in \ref Numerical_Flux_T.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Numerical_Flux_Input_T container.
 *  \return See brief. */
struct Numerical_Flux_Input_T* constructor_Numerical_Flux_Input_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Numerical_Flux_Input_T container.
void destructor_Numerical_Flux_Input_T
	(struct Numerical_Flux_Input_T* num_flux_i ///< Standard.
	);

/** \brief Constructor for a \ref Numerical_Flux_T container.
 *  \return See brief. */
struct Numerical_Flux_T* constructor_Numerical_Flux_T
	(const struct Numerical_Flux_Input_T* num_flux_i ///< \ref Numerical_Flux_Input_T.
	);

/// \brief Destructor for a \ref Numerical_Flux_T container.
void destructor_Numerical_Flux_T
	(struct Numerical_Flux_T* num_flux ///< Standard.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr_T for 1st order numerical fluxes only.
void compute_Numerical_Flux_1_T
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr_T for 1st and 2nd order numerical fluxes.
void compute_Numerical_Flux_12_T
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);
