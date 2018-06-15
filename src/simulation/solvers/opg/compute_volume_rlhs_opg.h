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

#ifndef DPG__compute_volume_rlhs_opg_h__INCLUDED
#define DPG__compute_volume_rlhs_opg_h__INCLUDED
/** \file
 *  \brief Provides real functions used for computing the volume contributions to the right and left-hand side
 *         (rlhs) terms of the OPG scheme.
 */

#include "def_templates_type_d.h"
#include "def_templates_compute_volume_rlhs_opg.h"
#include "compute_volume_rlhs_opg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_compute_volume_rlhs_opg.h"

#include <stdbool.h>

struct Flux_Ref;
struct Solver_Volume;

/** \brief Update the values of \ref Solver_Volume_T::sol_coef based on the updated \ref Solver_Volume_T::test_s_coef
 *         values. */
void update_coef_s_v_opg
	(const struct Simulation*const sim ///< Standard.
	);

/** \brief Constructor for the 'test' function 'diff'erentiation 'op'erator for the '1'st order 'v'olume term for the
 *         OPG scheme.
 *  \return See brief. */
const struct const_Matrix_d* constructor_test_diff_op_1v_opg
	(const struct Flux_Ref*const flux_r,     ///< Standard.
	 const struct Solver_Volume*const s_vol, ///< Standard.
	 const bool include_det_j                /**< Flag for whether the inverse Jacobian determinant should be
	                                          *   included. */
	);

#endif // DPG__compute_volume_rlhs_opg_h__INCLUDED
