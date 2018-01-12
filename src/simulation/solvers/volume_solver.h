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

#ifndef DPG__volume_solver_h__INCLUDED
#define DPG__volume_solver_h__INCLUDED
/** \file
 *  \brief Provides the interface for the real \ref Solver_Volume_T container and associated functions.
 */

#include "volume.h"
#include "geometry_blended.h"

#include "def_templates_type_d.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"
#include "def_templates_geometry.h"
#include "def_templates_volume_solver.h"
#include "volume_solver_T.h"
#include "undef_templates_type.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
#include "undef_templates_geometry.h"
#include "undef_templates_volume_solver.h"

struct Solver_Volume_c;

/// \brief Copy members from a real to a complex \ref Solver_Volume_T.
void copy_members_r_to_c_Solver_Volume
	(struct Solver_Volume_c*const s_vol,       ///< The complex \ref Solver_Volume_T.
	 const struct Solver_Volume*const s_vol_r, ///< The real \ref Solver_Volume_T.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

#endif // DPG__volume_solver_h__INCLUDED
