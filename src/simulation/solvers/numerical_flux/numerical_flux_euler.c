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

#include "multiarray.h"

#include "const_cast.h"
#include "flux.h"
#include "flux_euler.h"
#include "numerical_flux.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"

#include "def_templates_boundary.h"
#include "def_templates_const_cast_d.h"
#include "def_templates_flux.h"
#include "def_templates_math_d.h"
#include "def_templates_multiarray_d.h"
#include "def_templates_numerical_flux.h"

#include "numerical_flux_euler_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
