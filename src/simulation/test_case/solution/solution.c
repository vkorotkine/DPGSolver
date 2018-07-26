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

#include "solution.h"

#include "definitions_test_case.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "face_solver.h"
#include "element.h"
#include "element_solution.h"
#include "element_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "compute_volume_rlhs.h"
#include "flux.h"
#include "geometry.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solution_T.c"
#include "undef_templates_type.h"

struct Multiarray_d* constructor_rhs_v (const struct Simulation* sim, struct Solver_Volume* s_vol, const char node_kind)
{
	assert(list_is_derived_from("solver",'e',sim));
	assert((node_kind == 'c')); // Add support for other node kinds if required.

	assert(s_vol->rhs_0 != NULL);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	// If there is some odd behaviour when using this function, think further about whether the rhs term is a
	// coefficient or a value. Based on the case of identity mass matrix for the explicit scheme, it seems that rhs is
	// a coefficient.
	const struct Operator*const cv0_vt_vc = get_operator__cv0_vt_vc(s_vol);

	const struct const_Multiarray_d*const rhs_coef = (struct const_Multiarray_d*) s_vol->rhs_0;

	const ptrdiff_t ext_0 = cv0_vt_vc->op_std->ext_0,
	                ext_1 = rhs_coef->extents[1];

	struct Multiarray_d* rhs_v = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,ext_1}); // returned
	mm_NN1C_Operator_Multiarray_d(cv0_vt_vc,rhs_coef,rhs_v,op_format,rhs_coef->order,NULL,NULL);

	return rhs_v;
}

struct Multiarray_d* constructor_s_coef_bezier
	(const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	UNUSED(sim);
	const char op_format = get_set_op_format(0);

	const struct Volume*const vol         = (struct Volume*) s_vol;
	const struct Solver_Element*const s_e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref;
	const struct Operator*const ccSB0_vs_vs = get_Multiarray_Operator(s_e->ccSB0_vs_vs,(ptrdiff_t[]){0,0,p,p});

	struct Multiarray_d*const s_coef = s_vol->sol_coef;
	return constructor_mm_NN1_Operator_Multiarray_d(ccSB0_vs_vs,s_coef,'C',op_format,s_coef->order,NULL);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
