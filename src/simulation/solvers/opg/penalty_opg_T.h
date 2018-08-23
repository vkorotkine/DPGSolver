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
 *  \brief Provides functions used for the penalty term for both right and left-hand side (rlhs) terms of the OPG scheme.
 */

#include "def_templates_penalty_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver.h"

struct Flux_Ref_T;
struct Numerical_Flux_T;
struct Solver_Face_T;
struct Solver_Storage_Implicit;
struct Intrusive_List;

/** \brief Reset the values of \ref OPG_Solver_Face_T::bc_test_s and \ref OPG_Solver_Volume_T::bc_test_s for the current
 *         solver step.
 *
 *  Resetting these values is required as they are modified in the functions applying test boundary conditions such that
 *  the problem is neither over or under constrained.
 */
void reset_penalty_indicators_opg_T
	(const struct Intrusive_List*const faces ///< The list of faces for which to reset parameters.
	 );

/** \brief Version of \ref compute_rlhs_opg_f_fptr_T hitting a EXIT_UNSUPPORTED error.
 *
 *  This function is used for \ref OPG_Solver_Face_T::constructor_rlhs_penalty for interior faces.
 */
void constructor_rlhs_f_test_penalty_unsupported_T
	(const struct Flux_Ref_T*const flux_r,         ///< See brief.
	 const struct Numerical_Flux_T*const num_flux, ///< See brief.
	 struct Solver_Face_T*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi      ///< See brief.
	 );

/// \brief Version of \ref compute_rlhs_opg_f_fptr_T for rhs terms which does nothing.
void constructor_rlhs_f_test_penalty_do_nothing_T
	(const struct Flux_Ref_T*const flux_r,         ///< See brief.
	 const struct Numerical_Flux_T*const num_flux, ///< See brief.
	 struct Solver_Face_T*const s_face,            ///< See brief.
	 struct Solver_Storage_Implicit*const ssi      ///< See brief.
	 );

#include "undef_templates_penalty_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver.h"
