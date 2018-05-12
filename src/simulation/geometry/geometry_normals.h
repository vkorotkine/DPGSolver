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

#ifndef DPG__geometry_normals_h__INCLUDED
#define DPG__geometry_normals_h__INCLUDED
/** \file
 *  \brief Provides the interface to real functions from \ref geometry_normals_T.h.
 */

#include "def_templates_type_d.h"
#include "def_templates_geometry_normals.h"
#include "def_templates_face_solver.h"
#include "def_templates_multiarray.h"
#include "geometry_normals_T.h"
#include "undef_templates_type.h"
#include "undef_templates_geometry_normals.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_multiarray.h"

#include <stdbool.h>

/** \brief Return whether exact normals should be used.
 *  \return See brief. */
bool using_exact_normals ( );

#endif // DPG__geometry_normals_h__INCLUDED
