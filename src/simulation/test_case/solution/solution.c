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

#include "multiarray.h"
#include "matrix.h"

#include "computational_elements.h"
#include "face_solver.h"
#include "element.h"
#include "element_solution.h"
#include "element_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "flux.h"
#include "geometry.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_T.c"

struct Multiarray_d* constructor_rhs_v (const struct Simulation* sim, struct Solver_Volume* s_vol, const char node_kind)
{
	assert(list_is_derived_from("solver",'e',sim));
	assert((node_kind == 'c')); // Add support for other node kinds if required.

	assert(s_vol->rhs != NULL);

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int curved = vol->curved,
	          p      = s_vol->p_ref;

	// If there is some odd behaviour when using this function, think further about whether the rhs term is a
	// coefficient or a value. Based on the case of identity mass matrix for the explicit scheme, it seems that rhs is
	// a coefficient.
	const struct Operator* cv0_vs_vX = get_Multiarray_Operator(s_e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const struct const_Multiarray_T*const rhs_coef = (struct const_Multiarray_T*) s_vol->rhs;

	const ptrdiff_t ext_0 = cv0_vs_vX->op_std->ext_0,
	                ext_1 = rhs_coef->extents[1];

	struct Multiarray_T* rhs_v = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){ext_0,ext_1}); // returned
	mm_NN1C_Operator_Multiarray_T(cv0_vs_vX,rhs_coef,rhs_v,op_format,rhs_coef->order,NULL,NULL);

	return rhs_v;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
