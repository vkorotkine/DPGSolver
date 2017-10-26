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
 *  \brief Provides containers and functions relating to fluxes of the supported PDEs.
 */

struct Numerical_Flux_Input;
struct mutable_Numerical_Flux;
struct Simulation;
struct Face;

#include <stdbool.h>
#include "numerical_flux_euler.h"

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

/** \brief Function pointer to functions constructing the solution/gradients needed by the numerical flux at the face
 *         cubature nodes.
 *  \return Standard.
 *
 *  \param face \ref Face.
 *  \param sim  \ref Simulation.
 */
typedef const struct const_Multiarray_d* (*constructor_sg_fc_fptr)
	(const struct Face* face,
	 const struct Simulation* sim
	);

/// \brief Container holding data used for computing the numerical fluxes and their Jacobians.
struct Numerical_Flux_Input {
	const bool* compute_member; ///< Array of flags for which of the \ref Numerical_Flux members should be computed.

	const int d,     ///< \ref Simulation::d.
	          n_eq,  ///< \ref Test_Case::n_eq.
	          n_var; ///< \ref Test_Case::n_var.

	const bool has_1st_order, ///< \ref Test_Case::has_1st_order.
	           has_2nd_order; ///< \ref Test_Case::has_2nd_order.

	const struct const_Multiarray_d* normals_l; ///< The unit normal vector components as seen from the left.
	const struct const_Multiarray_d* xyz_l;     ///< The xyz coordinates as seen from the left.

	/** \brief Container for information from either side of the \ref Face.
	 *
	 *  The information for the first index `neigh_info[0]` relates to the what is termed the left side (whose outward
	 *  pointing unit normal is provided).
	 */
	struct Neigh_Info_NFI {
		const struct const_Multiarray_d* s; ///< The solution variables.
		const struct const_Multiarray_d* g; ///< The solution gradient variables.
	} neigh_info[2]; ///< \ref Neigh_Info_NFI.

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

/** \brief Version of \ref constructor_sg_fc_fptr constructing the solution needed by the numerical flux on the left
 *         side of the face at the face cubature nodes as seen from the left volume.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_s_l_fcl_interp
	(const struct Face* face,     ///< Defined for \ref constructor_sg_fc_fptr.
	 const struct Simulation* sim ///< Defined for \ref constructor_sg_fc_fptr.
	);

/** \brief Version of \ref constructor_sg_fc_fptr constructing the solution needed by the numerical flux on the right
 *         side of the face at the face cubature nodes as seen from the left volume.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_s_r_fcl_interp
	(const struct Face* face,     ///< Defined for \ref constructor_sg_fc_fptr.
	 const struct Simulation* sim ///< Defined for \ref constructor_sg_fc_fptr.
	);

/** \brief Version of \ref constructor_sg_fc_fptr returning NULL.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_sg_fc_null
	(const struct Face* face,     ///< Defined for \ref constructor_sg_fc_fptr.
	 const struct Simulation* sim ///< Defined for \ref constructor_sg_fc_fptr.
	);

/// \brief Destructor for a \ref Numerical_Flux container.
void destructor_Numerical_Flux
	(struct Numerical_Flux* num_flux ///< Standard.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr for 1st order numerical fluxes only.
void compute_Numerical_Flux_1
	(const struct Numerical_Flux_Input* num_flux_i, ///< Defined for \ref compute_Numerical_Flux_fptr.
	 struct mutable_Numerical_Flux* num_flux        ///< Defined for \ref compute_Numerical_Flux_fptr.
	);

/// \brief Version of \ref compute_Numerical_Flux_fptr for 1st and 2nd order numerical fluxes.
void compute_Numerical_Flux_12
	(const struct Numerical_Flux_Input* num_flux_i, ///< Defined for \ref compute_Numerical_Flux_fptr.
	 struct mutable_Numerical_Flux* num_flux        ///< Defined for \ref compute_Numerical_Flux_fptr.
	);

#endif // DPG__numerical_flux_h__INCLUDED
