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

#include "compute_volume_rlhs_fr_split_form.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "face_solver_fr_split_form.h"
#include "volume_solver_fr_split_form.h"
#include "element_solver_fr_split_form.h"

#include "compute_rlhs.h"
#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_fr_split_form.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the rhs and lhs terms for 1st order equations only.
/* static void compute_rlhs_1 */
/* 	(const struct Flux_Ref*const flux_r,      ///< See brief. */
/* 	 struct Solver_Volume*const s_vol,        ///< See brief. */
/* 	 struct Solver_Storage_Implicit*const ssi ///< See brief. */
/* 	); */



// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_volume_rlhs_fr_split_form_T.c"
#include "undef_templates_type.h"

/* #include "def_templates_type_dc.h" */
/* #include "compute_volume_rlhs_fr_split_form_T.c" */
/* #include "undef_templates_type.h" */

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_rlhs_v_fptr_T computing the lhs term for linearization wrt the solution.



// Level 1 ********************************************************************************************************** //

