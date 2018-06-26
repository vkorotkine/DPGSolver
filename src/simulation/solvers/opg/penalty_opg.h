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

#ifndef DPG__penalty_opg_h__INCLUDED
#define DPG__penalty_opg_h__INCLUDED
/** \file
 *  \brief Provides instantiations for \ref penalty_opg_T.h.
 */

#include "penalty_opg_advection.h"

#include "def_templates_type_d.h"
#include "penalty_opg_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "penalty_opg_T.h"
#include "undef_templates_type.h"

#define PENALTY_SCALING_OPG 1e4 ///< The scaling for the penalty term for OPG scheme test function boundary conditions.

#endif // DPG__penalty_opg_h__INCLUDED
