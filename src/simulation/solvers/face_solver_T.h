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
 *  \brief Provides the interface for the templated \ref Solver_Face_T container and associated functions.
 */

#include <stddef.h>

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_boundary.h"
#include "def_templates_face_solver.h"

/** \brief Container for data relating to the solver faces.
 *  \note Certain members are declared `const` despite requiring modification for adaptive simulations. Only members
 *        changing with every solver iteration are `mutable`.
 */
struct Solver_Face_T {
	struct Face face; ///< The base \ref Face.

	/// The index of the first degree of freedom (dof) of the face in relation to the global dof.
	const ptrdiff_t ind_dof;

	/// The reference order of the face. Need not be equal to the order of the solution in the face.
	const int p_ref;

	const int ml; ///< The mesh level of the face.

	/** Type of cubature to be used for the face. Options: 's'traight, 'c'urved. Set to curved whenever an adjacent
	 *  volume is curved. */
	const char cub_type;

	struct Multiarray_T* nf_coef; ///< The coefficients of the normal flux in the \ref Simulation::basis_sol.
	struct Multiarray_T* s_coef;  ///< The coefficients of the solution in the \ref Simulation::basis_sol.

	/// Values of the physical xyz coordinates at the face cubature nodes.
	const struct const_Multiarray_T* xyz_fc;

	/// Values of the physical xyz coordinates at the face cubature nodes on the exact domain boundary.
	const struct const_Multiarray_T* xyz_fc_ex_b;

	/// Values of the outward pointing unit normal vector at the face cubature nodes.
	const struct const_Multiarray_T* normals_fc;

	/** Values of the exact outward pointing unit normal vector at the face cubature nodes.
	 *
	 *  Unused unless one of the exact normal flags is enabled (boundary, all).
	 */
	const struct const_Multiarray_T* normals_fc_exact;

	/** The determinant of the face geometry Jacobian evaluated at the face cubature nodes. See (eq. (B.6),
	 *  \cite Zwanenburg2016) for the precise definition. */
	const struct const_Multiarray_T* jacobian_det_fc;

/// \todo Delete \ref Solver_Face_T::vol_jacobian_det_fc and \ref Solver_Face_T::metrics_fc if unused.
	/** The determinant of the volume geometry Jacobian evaluated at the face cubature nodes. This is the same as \ref
	 *  Solver_Volume_T::jacobian_det_vc but evaluated at the face cubature nodes. */
	const struct const_Multiarray_T* vol_jacobian_det_fc;

	/// The volume metric terms interpolated to the face cubature nodes.
	const struct const_Multiarray_T* metrics_fc;

	/// As for \ref Solver_Face_T::normals_fc but corresponding to the p1 geometry.
	const struct const_Multiarray_T* normals_p1;

	/// As for \ref Solver_Face_T::jacobian_det_fc but corresponding to the p1 geometry.
	const struct const_Multiarray_T* jacobian_det_p1;

	/** Construct 'r'ight numerical flux input members at face cubature nodes as seen from the left volume.
	 *
	 *  In the case of using the DG scheme, this function only computes the weak gradient contribution when the face
	 *  is on a domain boundary. This is necessary as the gradient must be set as the partially corrected value,
	 *  requiring \ref DG_Solver_Volume_T::grad_coef_v and DG_Solver_Face_T::Neigh_Info_DG::grad_coef_f, for internal
	 *  faces.
	 */
	constructor_Boundary_Value_fptr_T constructor_Boundary_Value_fcl;


	// Used exclusively for testing purposes.

	/// Exact values of the normal flux on the boundary at the face cubature nodes if the solution is known.
	const struct const_Multiarray_T* nf_fc;
};

/// \brief Constructor for a derived \ref Solver_Face_T.
void constructor_derived_Solver_Face_T
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref Solver_Face_T.
void destructor_derived_Solver_Face_T
	(struct Face* face_ptr ///< Pointer to the face.
	);

/** \brief Set the function pointers to the appropriate functions to compute values needed for the numerical flux
 *         computation. */
void set_function_pointers_face_num_flux_T
	(struct Solver_Face_T* s_face,  ///< Pointer to the \ref Solver_Face_T.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Get the appropriate sub-range of the \ref Solver_Element::w_fc operators.
 *  \return See brief. */
const struct const_Vector_R* get_operator__w_fc__s_e_T
	(const struct Solver_Face_T*const s_face ///< The current face.
	);

/** \brief Constructor for the mass matrix of the input face.
 *  \return See brief. */
const struct const_Matrix_T* constructor_mass_face_T
	(const struct Solver_Face_T*const s_face ///< Standard.
	 );

/** \brief Get the pointer to the appropriate \ref OPG_Solver_Element::cv0_vt_fc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vt_fc_T
	(const int side_index,                   ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face ///< Standard.
		);

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_boundary.h"
#include "undef_templates_face_solver.h"
