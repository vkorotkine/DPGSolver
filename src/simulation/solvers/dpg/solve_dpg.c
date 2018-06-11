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

#include "solve_dpg.h"

#include <assert.h>

#include "macros.h"

#include "face_solver.h"
#include "volume_solver.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_all_rlhs_dpg.h"
#include "compute_source_rlhs_dg.h"
#include "const_cast.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_dpg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solve_dpg_T.c"
#include "undef_templates_type.h"

double compute_rlhs_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	compute_all_rlhs_dpg(sim,ssi,sim->volumes);
	return compute_max_rhs_from_ssi(ssi);
}

void compute_flux_imbalances_dpg (struct Simulation*const sim)
{
	compute_flux_imbalances_faces_dpg(sim);
	compute_flux_imbalances_source_dg(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
