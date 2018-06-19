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
 *  \brief Provides templated functions used for computing the volume contributions to the right and left-hand side
 *         (rlhs) terms of the OPG scheme.
 */

#include "def_templates_compute_volume_rlhs_opg.h"
#include "def_templates_compute_rlhs.h"
#include "def_templates_matrix.h"
#include "def_templates_volume_solver_opg.h"

#include <stdbool.h>

struct Flux_Ref_T;
struct Solver_Volume_T;
struct OPG_Solver_Volume_T;
struct Simulation;
struct Solver_Storage_Implicit;
struct Intrusive_List;

/// \brief Compute the volume contributions to the rhs and lhs terms for the OPG scheme.
void compute_volume_rlhs_opg_T
	(const struct Simulation* sim,        ///< \ref Simulation.
	 struct Solver_Storage_Implicit* ssi, ///< \ref Solver_Storage_Implicit.
	 struct Intrusive_List* volumes       ///< The list of volumes.
	);

/** \brief Constructor for the operator used to compute \ref Solver_Volume_T::sol_coef from
 *         \ref Solver_Volume_T::test_s_coef for the OPG scheme.
 *  \return See brief. */
const struct const_Matrix_T* constructor_operator__test_s_coef_to_sol_coef_T
	(const struct Flux_Ref_T*const flux_r,             ///< Standard.
	 const struct OPG_Solver_Volume_T*const opg_s_vol, ///< Standard.
	 const bool include_m_inv                          /**< Flag for whether the inverse mass matrix should be
	                                                    *   included. */
		);

/** \brief Update the values of \ref Solver_Volume_T::sol_coef based on the updated \ref Solver_Volume_T::test_s_coef
 *         values. */
void update_coef_s_v_opg_T
	(const struct Simulation*const sim,  ///< Standard.
	 struct Intrusive_List*const volumes ///< The list of volumes for which to update the coefficients.
	 );

#include "undef_templates_compute_volume_rlhs_opg.h"
#include "undef_templates_compute_rlhs.h"
#include "undef_templates_matrix.h"
#include "undef_templates_volume_solver_opg.h"
