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

#include "compute_volume_rlhs_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "volume_solver_dg.h"
#include "element_solver_dg.h"

#include "compute_volume_rlhs.h"
#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref compute_rlhs_dg_fptr_T computing the rhs and lhs terms for any combination of 1st/2nd order
 *         equations only. */
static void compute_rlhs
	(const struct Flux_Ref* flux_r,       ///< See brief.
	 struct DG_Solver_Volume* dg_s_vol,   ///< See brief.
	 struct Solver_Storage_Implicit* ssi, ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "compute_volume_rlhs_dg_T.c"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void compute_rlhs
	(const struct Flux_Ref* flux_r, struct DG_Solver_Volume* dg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
/// \todo Add special case for collocated.
// If collocation is enabled, note that the diagonal weight scaling must be added back in to recover the symmetry of the
// residual Jacobian. Add it just before adding the contribution to the petsc mat. Also add for face terms and RHS
// terms (volume, face, source or simply the complete rhs).
assert(sim->collocated == false); // Add support in future.
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Solver_Volume* s_vol = (struct Solver_Volume*) dg_s_vol;
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc(s_vol);

	// rhs
	for (ptrdiff_t dim = 0; dim < DIM; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vt_vc.data[dim],flux_r->fr,dg_s_vol->rhs,op_format,2,&dim,NULL);
//print_Multiarray_d(dg_s_vol->rhs);

	// lhs
	struct Matrix_d* lhs = constructor_lhs_v_1(flux_r,s_vol,sim); // destructed
	set_petsc_Mat_row_col(ssi,s_vol,0,s_vol,0);
	add_to_petsc_Mat(ssi,(struct const_Matrix_d*)lhs);
	destructor_Matrix_d(lhs);
}
