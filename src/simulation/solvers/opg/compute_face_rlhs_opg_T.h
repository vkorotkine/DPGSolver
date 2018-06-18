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
 *  \brief Provides templated functions used for computing the face contributions to the right and left-hand side (rlhs)
 *         terms of the OPG scheme.
 *
 *  Please see the comments in \ref compute_face_rlhs_dg_T.h for the notation convention for the cross-terms in the
 *  linearized equations.
 */

#include "def_templates_compute_face_rlhs_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"
#include "def_templates_face_solver_opg.h"
#include "def_templates_matrix.h"
#include "def_templates_vector.h"

struct Simulation;
struct Solver_Storage_Implicit;
struct Intrusive_List;
struct Solver_Face_T;
struct OPG_Solver_Face_T;
struct Numerical_Flux_T;
struct Flux_Ref_T;

/** \brief Pointer to the function used to evaluate the rhs (and optionally lhs) face terms allowing for both reference
 *         flux and numerical flux inputs.
 *  \note The \ref Flux_Ref_T input is required for the linearization of the solution with respect to the test function.
 *
 *  \param flux_r   \ref Flux_Ref_T.
 *  \param num_flux \ref Numerical_Flux_T.
 *  \param s_face   \ref Solver_Face_T.
 *  \param ssi      \ref Solver_Storage_Implicit.
 */
typedef void (*compute_rlhs_opg_f_fptr_T)
	(const struct Flux_Ref_T*const flux_r,
	 const struct Numerical_Flux_T*const num_flux,
	 struct Solver_Face_T*const s_face,
	 struct Solver_Storage_Implicit*const ssi
		);

/// \brief Compute the face contributions to the rhs (and optionally lhs) terms for the DG scheme.
void compute_face_rlhs_opg_T
	(const struct Simulation*const sim,        ///< Standard.
	 struct Solver_Storage_Implicit*const ssi, ///< Standard.
	 struct Intrusive_List*const faces         ///< The list of faces.
	);

/** \brief Update the values of \ref Solver_Face_T::nf_coef based on the updated \ref Solver_Volume_T::test_s_coef
 *         values. */
void update_coef_nf_f_opg_T
	(const struct Simulation*const sim ///< Standard.
		);

/// \brief Container for operators needed for the assembly of LHS terms for the OPG scheme.
struct Lhs_Operators_OPG_T {
	const struct const_Vector_T* wJ_fc; ///< Face cubature weights "dot-multiplied" by the Jacobian determinant.

	/** 'c'oefficient to 'v'alue operators from the 'v'olume 't'est basis to 'f'ace 'c'ubature nodes.
	 *
	 *  The indices are used to denote the following operators:
	 *  - [0]: ll operator ('l'eft  basis -> 'l'eft  nodes);
	 *  - [1]: rl operator ('r'ight basis -> 'l'eft  nodes; permutation of node ordering);
	 */
	const struct const_Matrix_R* cv0_vt_fc[2];
};

/** \brief Constructor for a \ref Lhs_Operators_OPG container.
 *  \return See brief. */
const struct Lhs_Operators_OPG_T* constructor_Lhs_Operators_OPG_T
	(const struct OPG_Solver_Face_T*const opg_s_face ///< Standard.
	);

/// \brief Destructor for a \ref Lhs_Operators_OPG container.
void destructor_Lhs_Operators_OPG_T
	(const struct Solver_Face_T*const s_face, ///< Standard.
	 const struct Lhs_Operators_OPG_T* ops    ///< Standard.
	);

#include "undef_templates_compute_face_rlhs_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_face_solver_opg.h"
#include "undef_templates_matrix.h"
#include "undef_templates_vector.h"
