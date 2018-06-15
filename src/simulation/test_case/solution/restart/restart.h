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

#ifndef DPG__restart_h__INCLUDED
#define DPG__restart_h__INCLUDED
/** \file
 *  \brief Provides the real interface to \ref restart_T.h
 */

#include "def_templates_type_d.h"
#include "restart_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "restart_T.h"
#include "undef_templates_type.h"

/** \brief Return a stack allocated `char*` holding the name (including full path) of the restart file.
 *  \return See brief. */
const char* get_restart_name ( );

#endif // DPG__restart_h__INCLUDED
