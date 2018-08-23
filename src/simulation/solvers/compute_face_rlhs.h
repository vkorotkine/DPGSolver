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

#ifndef DPG__compute_face_rlhs_h__INCLUDED
#define DPG__compute_face_rlhs_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the face contributions to the right and left-hand side (rlhs) terms
 *         of supported schemes.
 */

#include "def_templates_type_d.h"
#include "compute_face_rlhs_T.h"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "compute_face_rlhs_T.h"
#include "undef_templates_type.h"

struct Numerical_Flux;
struct Solver_Face;

/// \brief Version of \ref scale_by_Jacobian_fptr_T for 1st order equations only.
void scale_by_Jacobian_i1
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

/// \brief Version of \ref scale_by_Jacobian_fptr_T for 2nd order equations only.
void scale_by_Jacobian_i2
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

/// \brief Version of \ref scale_by_Jacobian_fptr_T for both 1st and 2nd order equations.
void scale_by_Jacobian_i12
	(struct Numerical_Flux*const num_flux, ///< See brief.
	 const struct Solver_Face*const s_face ///< See brief.
		);

#endif // DPG__compute_face_rlhs_h__INCLUDED
