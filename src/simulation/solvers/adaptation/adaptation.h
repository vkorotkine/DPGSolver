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

#ifndef DPG__adaptation_h__INCLUDED
#define DPG__adaptation_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used to perform hp adaption.
 */

#include <stdbool.h>

struct Simulation;
struct Solver_Face;

/// \brief Container for data relating to element adaptation in the domain.
struct Adaptation_Data {
	int adapt_h[2]; ///< Number of uniform refinement [0] and coarsening [1] adaptations to perform.

	/// Coordinates of vertices for which adjacent elements should be refined.
	const struct const_Multiarray_d* xyz_ve_refine;

	const struct const_Vector_i* xyz_ve_ml; ///< Mesh levels to which xyz_ve_refine should be refined.
	const struct const_Vector_i* xyz_ve_p;  ///< Polynomial order to which xyz_ve_refine should be refined/coarsened.
};

/// \brief Adapt the computational elements using the input adaptation strategy.
void adapt_hp
	(struct Simulation* sim,                       ///< \ref Simulation.
	 const int adapt_strategy,                     /**< The adaptation strategy. Options:
	                                                *   see \ref definitions_adaptation.h. */
	 const struct Adaptation_Data*const adapt_data ///< \ref Adaptation_Data.
	);

/** \brief Constructor for the geometry values at the face geometry nodes interpolated from the volume of specified
 *         side_index including reordering if the destination side_index differs.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_geom_fg
	(const int side_index,                  ///< The side index of the neighbouring volume.
	 const int side_index_dest,             ///< The side index of the destination neighbouring volume.
	 const struct Solver_Face*const s_face, ///< The current \ref Solver_Face.
	 const bool use_pg_face,                /**< Flag for whether the face geometry order should be used instead of
	                                         *   the face reference polynomial order. */
	 const bool use_full_face               /**< Flag for whether the whole face (not the sub-face if applicable)
	                                         *   should be used. */
	);

#endif // DPG__adaptation_h__INCLUDED
