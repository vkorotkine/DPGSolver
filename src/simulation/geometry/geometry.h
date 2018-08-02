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

#ifndef DPG__geometry_h__INCLUDED
#define DPG__geometry_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for geometry processing.
 */

#include "def_templates_type_d.h"
#include "geometry_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "geometry_T.h"
#include "undef_templates_type.h"

#include "geometry_blended.h"
#include "geometry_parametric.h"
#include "geometry_surface.h"

/** \brief Return whether the internal faces of the mesh should use straight geometry as specified in the input file.
 *  \return See brief. */
bool is_internal_geom_straight ( );

/** \brief Return whether the geometry coefficients are dependent upon the pointers for faces adjacent to volumes having
 *         been previously set.
 *  \return See brief.
 *
 *  For example, for a parametric domain in which the internal faces are to be straight, the function used to correct
 *  the curved faces loops over faces adjacent to volumes and correcting the volume geometry if required.
 */
bool geometry_depends_on_face_pointers( );

#endif // DPG__geometry_h__INCLUDED
