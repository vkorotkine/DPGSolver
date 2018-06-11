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

#include "macros.h"

#include "def_templates_volume_solver_opg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"

#include "def_templates_volume_solver.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_OPG_Solver_Volume_T (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Solver_Volume_T*const s_vol         = (struct Solver_Volume_T*) volume_ptr;
	struct OPG_Solver_Volume_T*const opg_s_vol = (struct OPG_Solver_Volume_T*) volume_ptr;
	UNUSED(sim);

	opg_s_vol->m_inv = constructor_inverse_mass_T(s_vol,NULL); // destructed
}

void destructor_derived_OPG_Solver_Volume_T (struct Volume* volume_ptr)
{
	struct OPG_Solver_Volume_T* opg_s_vol = (struct OPG_Solver_Volume_T*) volume_ptr;

	destructor_const_Matrix_R(opg_s_vol->m_inv);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

#include "undef_templates_volume_solver_opg.h"

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"

#include "undef_templates_volume_solver.h"
