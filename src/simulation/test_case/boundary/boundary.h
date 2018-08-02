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

#ifndef DPG__boundary_h__INCLUDED
#define DPG__boundary_h__INCLUDED
/** \file
 *  \brief Provides containers and functions relating to boundary conditions of the supported PDEs.
 */

///\{ \name The maximum number of outputs from the boundary functions.
#define MAX_BV_OUT 6 ///< See the members of \ref Boundary_Value_T.
///\}

#include "def_templates_type_d.h"
#include "boundary_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "boundary_T.h"
#include "undef_templates_type.h"

/** \brief Return whether exact boundary normals should be used for boundary condition computation.
 *  \return See brief. */
bool using_exact_normals_for_boundary ( );

/** \brief Return `true` if the input boundary condition is adjoint consistent, `false` otherwise.
 *  \return See brief. */
bool using_adjoint_consistent_bc
	(const int bc ///< The input boundary condition.
		);

#endif // DPG__boundary_h__INCLUDED
