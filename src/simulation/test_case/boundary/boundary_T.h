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
 *  \brief Provides templated containers and functions relating to boundary conditions of the supported PDEs.
 */

#include <stdbool.h>

struct Boundary_Value_Input_T;
struct Boundary_Value_T;
struct Solver_Face_T;
struct Simulation;

/** \brief Function pointer to functions constructing the face-specific members \ref Boundary_Value_Input.
 *  \return Standard.
 *
 *  \param bv_i   \ref Boundary_Value_Input.
 *  \param s_face \ref Solver_Face.
 *  \param sim    \ref Simulation.
 */
typedef void (*constructor_Boundary_Value_Input_face_fptr_T)
	(struct Boundary_Value_Input_T* bv_i,
	 const struct Solver_Face_T* s_face,
	 const struct Simulation* sim
	);

/** \brief Function pointer to functions constructing the members \ref Boundary_Value either from the right volume or by
 *         calling a boundary condition function.
 *  \return Standard.
 *
 *  \param bv     \ref Boundary_Value.
 *  \param bv_i   \ref Boundary_Value_Input.
 *  \param s_face \ref Solver_Face.
 *  \param sim    \ref Simulation.
 */
typedef void (*constructor_Boundary_Value_fptr_T)
	(struct Boundary_Value_T* bv,
	 const struct Boundary_Value_Input_T* bv_i,
	 const struct Solver_Face_T* s_face,
	 const struct Simulation* sim
	);

/// \brief Container holding data used for computing the boundary condition values and their Jacobians.
struct Boundary_Value_Input_T {
	const char* input_path;     ///< Pointer to \ref Simulation::input_path.
	const bool* compute_member; ///< Flags for which of the \ref Boundary_Value members should be computed.

	const int n_eq,  ///< \ref Test_Case::n_eq.
	          n_var; ///< \ref Test_Case::n_var.

/// \todo geometry: Real -> Templated
	const struct const_Multiarray_R* normals; ///< The unit normal vector components.
	const struct const_Multiarray_R* xyz;     ///< The xyz coordinates.

	const struct const_Multiarray_T* s; ///< The solution variables.
	const struct const_Multiarray_T* g; ///< The solution gradient variables.
};

/// \brief Container storing the boundary condition values and their Jacobians.
struct Boundary_Value_T {
	const struct const_Multiarray_T* s; ///< The solution variables values on the boundary.
	const struct const_Multiarray_T* g; ///< The solution gradient variables values on the boundary.

	const struct const_Multiarray_T* ds_ds; ///< The Jacobian of the boundary solution wrt the solution.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_Boundary_Value_Input_face_fptr constructing only the solution using members from
 *         the face and interpolated from the left volume at the face cubature nodes as seen from the left volume.
 *  \return See brief. */
void constructor_Boundary_Value_Input_face_s_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, ///< Defined for \ref constructor_Boundary_Value_Input_face_fptr.
	 const struct Solver_Face_T* face,      ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Simulation* sim         ///< Defined for \ref constructor_Boundary_Value_Input_face_fptr.
	);

/// \brief Destructor for a \ref Boundary_Value_Input container.
void destructor_Boundary_Value_Input_T
	(struct Boundary_Value_Input_T* bv_i ///< Standard.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr interpolated from the right volume at the face cubature nodes
 *         as seen from the left volume.
 *  \return See brief. */
void constructor_Boundary_Value_s_fcl_interp_T
	(struct Boundary_Value_T* bv,               ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Boundary_Value_Input_T* bv_i, ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Solver_Face_T* face,            ///< Defined for \ref constructor_Boundary_Value_fptr.
	 const struct Simulation* sim               ///< Defined for \ref constructor_Boundary_Value_fptr.
	);

/// \brief Destructor for a \ref Boundary_Value container.
void destructor_Boundary_Value_T
	(struct Boundary_Value_T* bv ///< Standard.
	);
