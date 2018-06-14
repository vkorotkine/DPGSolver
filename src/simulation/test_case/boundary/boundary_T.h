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

#include "def_templates_multiarray.h"
#include "def_templates_boundary.h"
#include "def_templates_operators.h"
#include "def_templates_face_solver.h"

struct Boundary_Value_Input_T;
struct Boundary_Value_T;
struct Solver_Face_T;
struct Simulation;

/** \brief Function pointer to functions constructing the face-specific members \ref Boundary_Value_Input_T.
 *  \return Standard.
 *
 *  \param bv_i   \ref Boundary_Value_Input_T.
 *  \param s_face \ref Solver_Face_T.
 *  \param sim    \ref Simulation.
 */
typedef void (*constructor_Boundary_Value_Input_face_fptr_T)
	(struct Boundary_Value_Input_T* bv_i,
	 const struct Solver_Face_T* s_face,
	 const struct Simulation* sim
	);

/** \brief Function pointer to functions constructing the members \ref Boundary_Value_T either from the right volume or by
 *         calling a boundary condition function.
 *  \return Standard.
 *
 *  \param bv     \ref Boundary_Value_T.
 *  \param bv_i   \ref Boundary_Value_Input_T.
 *  \param s_face \ref Solver_Face_T.
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
	/** Flags for which of the \ref Boundary_Value_T members should be computed.
	 *  - [0]: s
	 *  - [1]: ds_ds
	 *  - [2]: g
	 *  - [3]: dg_dg
	 *  - [4]: dg_ds
	 *  - [5]: ds_dg
	 */
	const bool* compute_member;

	const int n_eq,  ///< \ref Test_Case_T::n_eq.
	          n_var; ///< \ref Test_Case_T::n_var.

	int bc;   ///< \ref Face::bc.
	double h; ///< \ref Face::h.
	int p;    ///< \ref Solver_Face_T::p_ref.

	const struct const_Multiarray_T* normals;     ///< The unit normal vector components.
	const struct const_Multiarray_T* normals_std; ///< Standard unit normal vector components (computed from metrics).
	const struct const_Multiarray_T* xyz;     ///< The xyz coordinates.
	const struct const_Multiarray_T* xyz_ex;  ///< \ref Solver_Face_T::xyz_fc_ex_b.

	/// \ref Solver_Face_T::jacobian_det_fc. \todo Remove if unused (was possibly used only for testing purposes).
	const struct const_Multiarray_T* jacobian_det_fc;

	const struct const_Multiarray_T* s; ///< The solution variables.
	const struct const_Multiarray_T* g; ///< The solution gradient variables.
};

/// \brief Container storing the boundary condition values and their Jacobians.
struct Boundary_Value_T {
	/** 'n'ormal 'f'lux for the 'E'nergy equation in cases where a constant value is used (e.g. Diabatic wall BCs).
	 *
	 *  In these cases, neither boundary values, gradients, nor their Jacobians are required as they are not used to
	 *  compute the normal numerical flux and they are all set to NULL. */
	Real nf_E;

	bool nf_E_provided; ///< Flag indicating whether \ref Boundary_Value_T::nf_E was set and should be used.

	const struct const_Multiarray_T* s; ///< The solution variables values on the boundary.
	const struct const_Multiarray_T* g; ///< The solution gradient variables values on the boundary.

	/** The Jacobian of the boundary solution wrt the solution. The ordering was chosen intuitively for consistency
	 *  with the non-linearized term and the treatment of other linearized terms throughout the code. Memory changes
	 *  fastest with respect to the boundary variables, then with respect to the derivative variables (i.e. the
	 *  columns of ds_ds are given by: drho_b/drho_i, drhou_b/drho_i, ..., drho_b/drhou_i, ... where "_b" and "_i"
	 *  represent 'b'oundary and 'i'nternal variables respectively). Note that "_b" and "_i" are generally replaced
	 *  with "_r" ('r'ight) and "_l" ('l'eft) throughout for consistency with naming conventions used for internal
	 *  faces.
	 */
	const struct const_Multiarray_T* ds_ds;

	/** The Jacobian of the boundary gradient wrt to the gradient. See comments for \ref Boundary_Value_T::ds_ds in
	 *  relation to the ordering (fastest on boundary variable/dimension, then on internal variable/dimension). */
	const struct const_Multiarray_T* dg_dg;

	/** The Jacobian of the boundary gradient wrt to the solution. See comments for \ref Boundary_Value_T::ds_ds in
	 *  relation to the ordering (fastest on boundary variable/dimension, then on internal variable). */
	const struct const_Multiarray_T* dg_ds;
};

/// \brief `mutable` version of \ref Boundary_Value_T.
struct mutable_Boundary_Value_T {
	Real nf_E;          ///< See brief.
	bool nf_E_provided; ///< See brief.

	struct Multiarray_T* s;     ///< See brief.
	struct Multiarray_T* g;     ///< See brief.
	struct Multiarray_T* ds_ds; ///< See brief.
	struct Multiarray_T* dg_dg; ///< See brief.
	struct Multiarray_T* dg_ds; ///< See brief.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_Boundary_Value_Input_face_fptr_T constructing only the solution using members from
 *         the face and interpolated from the left volume at the face cubature nodes as seen from the left volume.
 *  \return See brief. */
void constructor_Boundary_Value_Input_face_s_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* s_face,  ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/** \brief Version of \ref constructor_Boundary_Value_Input_face_fptr_T constructing the solution and solution gradient
 *         using members from the face and interpolated from the left volume at the face cubature nodes as seen from the
 *         left volume.
 *  \return See brief. */
void constructor_Boundary_Value_Input_face_sg_fcl_interp_T
	(struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* s_face,  ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/// \brief Destructor for a \ref Boundary_Value_Input_T container.
void destructor_Boundary_Value_Input_T
	(struct Boundary_Value_Input_T* bv_i ///< Standard.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T interpolated from the right volume at the face cubature nodes
 *         as seen from the left volume.
 *  \return See brief. */
void constructor_Boundary_Value_s_fcl_interp_T
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* s_face,        ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Destructor for a \ref Boundary_Value_T container.
void destructor_Boundary_Value_T
	(struct Boundary_Value_T* bv ///< Standard.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T constructing gradient components of the input \ref
 *         Boundary_Value_T in the case of the ghost values for the gradient being specified as the same as the internal
 *         values. */
void constructor_Boundary_Value_T_grad_from_internal
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const int n_var                            ///< The number of variables.
	);

/** \brief Constructor for the solution interpolating from the neighbouring volume to the face cubature nodes.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_s_fc_interp_T
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const face ///< Standard.
	);

/** \brief Constructor for the solution gradient interpolating from the neighbouring volume to the face cubature nodes.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_g_fc_interp_T
	(const int side_index,                 ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const face ///< Standard.
	);

#include "undef_templates_multiarray.h"
#include "undef_templates_boundary.h"
#include "undef_templates_operators.h"
#include "undef_templates_face_solver.h"
