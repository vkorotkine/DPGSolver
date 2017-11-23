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

#ifndef DPG__face_solver_h__INCLUDED
#define DPG__face_solver_h__INCLUDED
/** \file
 *  \brief Provides the interface for the \ref Solver_Face container and associated functions.
 */

#include <stddef.h>
#include "face.h"
#include "boundary.h"

/** \brief Container for data relating to the solver faces.
 *  \note Certain members are declared `const` despite requiring modification for adaptive simulations. Only members
 *        changing with every solver iteration are `mutable`.
 */
struct Solver_Face {
	struct Face face; ///< The base \ref Face.

	/// The index of the first degree of freedom (dof) of the face in relation to the global dof.
	const ptrdiff_t ind_dof;

	/// The reference order of the face. Need not be equal to the order of the solution in the face.
	const int p_ref;

	const int ml; ///< The mesh level of the face.

	/** Type of cubature to be used for the face. Options: 's'traight, 'c'urved. Set to curved whenever an adjacent
	 *  volume is curved. */
	const char cub_type;

	/// The coefficients of the normal flux in the \ref Simulation::basis_sol.
	struct Multiarray_d* nf_coef;

	/// Values of the physical xyz coordinates at the face cubature nodes.
	const struct const_Multiarray_d*const xyz_fc;

	/// Values of the outward pointing unit normal vector at the face cubature nodes.
	const struct const_Multiarray_d*const normals_fc;

	/** The determinant of the face geometry Jacobian evaluated at the face cubature nodes. See (eq. (B.6),
	 *  \cite Zwanenburg2016) for the precise definition. */
	const struct const_Multiarray_d*const jacobian_det_fc;

	/// Construct 'r'ight numerical flux input members at face cubature nodes as seen from the left volume.
	constructor_Boundary_Value_fptr constructor_Boundary_Value_fcl;
};

/// \brief Constructor for a derived \ref Solver_Face.
void constructor_derived_Solver_Face
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Solver_Face.
void destructor_derived_Solver_Face
	(struct Face* face_ptr ///< Pointer to the face.
	);

#endif // DPG__face_solver_h__INCLUDED
