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
 */

#include "volume_solver_dpg.h"

#include "element_solver_dpg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "compute_all_rlhs_dpg.h"
#include "compute_volume_rlhs.h"

#include "test_complex_volume_solver_dpg.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "volume_solver_dpg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "volume_solver_dpg_T.c"
#include "undef_templates_type.h"

void copy_members_r_to_c_DPG_Solver_Volume
	(struct DPG_Solver_Volume_c*const dpg_s_vol, const struct DPG_Solver_Volume*const dpg_s_vol_r,
	 const struct Simulation*const sim)
{
	UNUSED(sim);
	destructor_derived_DPG_Solver_Volume_T((struct Volume*)dpg_s_vol);
	dpg_s_vol->norm_op_H0 = constructor_copy_const_Matrix_d(dpg_s_vol_r->norm_op_H0); // destructed
	dpg_s_vol->norm_op_H1 = constructor_copy_const_Matrix_d(dpg_s_vol_r->norm_op_H1); // destructed
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
