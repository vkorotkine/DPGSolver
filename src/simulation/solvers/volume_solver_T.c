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

#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "definitions_mesh.h"


#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_geometry.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Set the appropriate xyz surface constructor geometry function pointer if applicable.
 *
 *  These functions serve to set the values of the geometry coordinates on curved element faces based on a given
 *  parametrization of the current geometry "patch". A patch boundary occurs whenever there exists a region of limited
 *  geometry smoothness (i.e. below C^{infinity}).
 */
static void set_function_pointers_constructor_xyz_surface
	(struct Solver_Volume_T*const s_vol, ///< The current volume.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Volume_T (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) volume_ptr;

	const_cast_ptrdiff(&s_vol->ind_dof,-1);
	const_cast_ptrdiff(&s_vol->ind_dof_constraint,-1);
	const_cast_i(&s_vol->p_ref,sim->p_ref[0]);
	const_cast_i(&s_vol->ml,0);

	const_constructor_move_Multiarray_R(&s_vol->geom_coef,constructor_default_Multiarray_R()); // destructed
	set_function_pointers_constructor_xyz_surface(s_vol,sim);

	s_vol->sol_coef  = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){0,0});   // destructed
	s_vol->grad_coef = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){0,0,0}); // destructed

	const_constructor_move_Multiarray_R(
		&s_vol->metrics_vm,constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_R(
		&s_vol->metrics_vc,constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){0,0,0}));  // destructed
	const_constructor_move_Multiarray_R(
		&s_vol->jacobian_det_vc,constructor_empty_Multiarray_R('C',1,(ptrdiff_t[]){0})); // destructed

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	s_vol->flux_imbalance = constructor_empty_Vector_T(test_case->n_var); // destructed
	s_vol->l_mult = constructor_zero_Multiarray_T('C',1,(ptrdiff_t[]){test_case->n_eq}); // destructed
	s_vol->rhs = NULL;
}

void destructor_derived_Solver_Volume_T (struct Volume* volume_ptr)
{
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) volume_ptr;

	destructor_const_Multiarray_R(s_vol->geom_coef);
	destructor_Multiarray_T(s_vol->sol_coef);
	destructor_Multiarray_T(s_vol->grad_coef);
	destructor_const_Multiarray_R(s_vol->metrics_vm);
	destructor_const_Multiarray_R(s_vol->metrics_vc);
	destructor_const_Multiarray_R(s_vol->jacobian_det_vc);
	destructor_Vector_T(s_vol->flux_imbalance);
	destructor_Multiarray_T(s_vol->l_mult);
	destructor_conditional_Multiarray_T(s_vol->rhs);
}

const struct const_Vector_d* get_operator__w_vc__s_e_T (const struct Solver_Volume_T* s_vol)
{
	struct Volume* vol               = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	return get_const_Multiarray_Vector_d(s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_constructor_xyz_surface
	(struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
{
	const struct Volume*const vol = (struct Volume*) s_vol;
	if (!vol->boundary || !vol->curved) {
		s_vol->constructor_xyz_surface = NULL;
		return;
	}

	if (sim->domain_type == DOM_PARAMETRIC) {
		s_vol->constructor_xyz_surface = set_constructor_xyz_surface_fptr_T(NULL,-1,sim->domain_type);
	} else if (sim->domain_type == DOM_BLENDED) {
		const struct Test_Case_T*const test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
		const int geom_par = test_case->geom_parametrization;
		if ((strcmp(sim->geom_name,"n-cylinder_hollow_section") == 0) ||
		    (strcmp(sim->geom_name,"n-cylinder_hollow")         == 0)) {
			s_vol->constructor_xyz_surface =
				set_constructor_xyz_surface_fptr_T("n-cylinder",geom_par,sim->domain_type);
		} else {
			EXIT_ERROR("Unsupported: %s, %s, %d\n",sim->geom_name,sim->geom_spec,sim->domain_type);
		}
	} else {
		EXIT_ERROR("Unsupported: %d\n",sim->domain_type);
	}
}
