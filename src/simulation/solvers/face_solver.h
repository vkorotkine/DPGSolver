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

#ifndef DPG__face_solver_h__INCLUDED
#define DPG__face_solver_h__INCLUDED
/** \file
 *  \brief Provides the interface for the real \ref Solver_Face_T container and associated functions.
 */

#include "face.h"
#include "boundary.h"

#include "def_templates_type_d.h"
#include "face_solver_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "face_solver_T.h"
#include "undef_templates_type.h"

struct Solver_Face_c;

/// \brief Copy members from a real to a complex \ref Solver_Face_T.
void copy_members_r_to_c_Solver_Face
	(struct Solver_Face_c*const s_face,       ///< The complex \ref Solver_Face_T.
	 const struct Solver_Face*const s_face_r, ///< The real \ref Solver_Face_T.
	 const struct Simulation*const sim        ///< \ref Simulation.
	);

/** \brief Return `true` if the face is conforming; `false` otherwise.
 *  \return See brief. */
bool face_is_conforming
	(const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

/** \brief Return the index of the side corresponding to the volume which is dominant in relation to the geometry
 *         specification.
 *  \return See brief.
 *
 *  The dominant volume is that which constrains the geometry the most.
 */
int get_dominant_geom_vol_side_index
	(const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

/** \brief Return the geometry polynomial order to be used for the face.
 *  \return See brief. */
int compute_face_geometry_order
	(const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

/** \brief Return the reference polynomial order for the face.
 *  \return See brief. */
int compute_face_reference_order
	(const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

#endif // DPG__face_solver_h__INCLUDED
