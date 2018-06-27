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

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_bc.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_boundary.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Set the function pointers to the appropriate functions to compute boundary values needed for the numerical
 *         flux computation. */
static void set_function_pointers_num_flux_bc
	(struct Solver_Face_T* s_face, ///< Defined for \ref set_function_pointers_face_num_flux_T.
	 const struct Simulation* sim  ///< Defined for \ref set_function_pointers_face_num_flux_T.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_Solver_Face_T (struct Face* face_ptr, const struct Simulation* sim)
{
	struct Solver_Face_T* s_face = (struct Solver_Face_T*) face_ptr;

	const_cast_ptrdiff(&s_face->ind_dof,-1);
	const_cast_i(&s_face->p_ref,sim->p_ref[0]);
	const_cast_i(&s_face->ml,0);
	const_cast_c(&s_face->cub_type,(check_for_curved_neigh((struct Face*)s_face) ? 'c' : 's'));

	s_face->nf_coef = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){0,0});   // destructed

	s_face->xyz_fc              = constructor_empty_const_Multiarray_T('C',2,(ptrdiff_t[]){0,0}); // destructed
	s_face->xyz_fc_ex_b         = constructor_empty_const_Multiarray_T('C',2,(ptrdiff_t[]){0,0}); // destructed
	s_face->normals_fc          = constructor_empty_const_Multiarray_T('R',2,(ptrdiff_t[]){0,0}); // destructed
	s_face->normals_fc_exact    = constructor_empty_const_Multiarray_T('R',2,(ptrdiff_t[]){0,0}); // destructed
	s_face->jacobian_det_fc     = constructor_empty_const_Multiarray_T('C',1,(ptrdiff_t[]){0});   // destructed
	s_face->vol_jacobian_det_fc = constructor_empty_const_Multiarray_T('C',1,(ptrdiff_t[]){0});   // destructed
	s_face->metrics_fc          = constructor_empty_const_Multiarray_T('C',3,(ptrdiff_t[]){0,0,0}); // destructed
	s_face->normals_p1          = constructor_empty_const_Multiarray_T('R',2,(ptrdiff_t[]){0,0}); // destructed
	s_face->jacobian_det_p1     = constructor_empty_const_Multiarray_T('C',1,(ptrdiff_t[]){0});   // destructed

	set_function_pointers_face_num_flux_T(s_face,sim);

	s_face->nf_fc = NULL;
}

void destructor_derived_Solver_Face_T (struct Face* face_ptr)
{
	struct Solver_Face_T* face = (struct Solver_Face_T*) face_ptr;

	destructor_Multiarray_T(face->nf_coef);

	destructor_const_Multiarray_T(face->xyz_fc);
	destructor_const_Multiarray_T(face->xyz_fc_ex_b);
	destructor_const_Multiarray_T(face->normals_fc);
	destructor_const_Multiarray_T(face->normals_fc_exact);
	destructor_const_Multiarray_T(face->jacobian_det_fc);
	destructor_const_Multiarray_T(face->vol_jacobian_det_fc);
	destructor_const_Multiarray_T(face->metrics_fc);
	destructor_const_Multiarray_T(face->normals_p1);
	destructor_const_Multiarray_T(face->jacobian_det_p1);

	destructor_conditional_const_Multiarray_T(face->nf_fc);
}

void set_function_pointers_face_num_flux_T (struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	const struct Face* face = (struct Face*) s_face;
	if (!face->boundary) {
		struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
		switch (test_case->pde_index) {
		case PDE_ADVECTION: // fallthrough
		case PDE_EULER:     // fallthrough
		case PDE_BURGERS_INVISCID:
			s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_s_fcl_interp_T;
			break;
		case PDE_DIFFUSION:     // fallthrough
		case PDE_NAVIER_STOKES:
			/// Note: Does not construct the solution gradient term. See comment in \ref Solver_Face_T.
			s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_s_fcl_interp_T;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
			break;
		}
	} else {
		set_function_pointers_num_flux_bc(s_face,sim);
	}
}

const struct const_Vector_R* get_operator__w_fc__s_e_T (const struct Solver_Face_T*const s_face)
{
	const int side_index = 0;
	const struct Face*const face          = (struct Face*) s_face;
	const struct Volume*const vol         = face->neigh_info[side_index].volume;
	const struct Solver_Element*const s_e = (struct Solver_Element*) vol->element;

	const int p_f = s_face->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	return get_const_Multiarray_Vector_d(s_e->w_fc[curved],(ptrdiff_t[]){0,0,0,0,p_f,p_f});
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref set_function_pointers_num_flux_bc for the linear advection equation.
static void set_function_pointers_num_flux_bc_advection
	(struct Solver_Face_T* s_face ///< See brief.
	);

/// \brief Version of \ref set_function_pointers_num_flux_bc for the diffusion equation.
static void set_function_pointers_num_flux_bc_diffusion
	(struct Solver_Face_T* s_face ///< See brief.
	);

/// \brief Version of \ref set_function_pointers_num_flux_bc for the Euler equations.
static void set_function_pointers_num_flux_bc_euler
	(struct Solver_Face_T* s_face ///< See brief.
	);

/// \brief Version of \ref set_function_pointers_num_flux_bc for the Navier-Stokes equations.
static void set_function_pointers_num_flux_bc_navier_stokes
	(struct Solver_Face_T* s_face ///< See brief.
	);

static void set_function_pointers_num_flux_bc (struct Solver_Face_T* s_face, const struct Simulation* sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	switch (test_case->pde_index) {
	case PDE_ADVECTION:     set_function_pointers_num_flux_bc_advection(s_face); break;
	case PDE_DIFFUSION:     set_function_pointers_num_flux_bc_diffusion(s_face); break;
	case PDE_EULER:         set_function_pointers_num_flux_bc_euler(s_face);     break;
	case PDE_NAVIER_STOKES: set_function_pointers_num_flux_bc_navier_stokes(s_face); break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
}

// Level 1 ********************************************************************************************************** //

static void set_function_pointers_num_flux_bc_advection (struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_INFLOW: case BC_INFLOW_ALT1: case BC_INFLOW_ALT2:
	case BC_OUTFLOW: case BC_OUTFLOW_ALT1: case BC_OUTFLOW_ALT2:
	case BC_UPWIND: case BC_UPWIND_ALT1: case BC_UPWIND_ALT2:
	case BC_UPWIND_ALT3: case BC_UPWIND_ALT4: case BC_UPWIND_ALT5:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_advection_upwind;
		break;
	case BC_SLIPWALL:
		printf("This BC was deprecated and removed from the code.\n");
		printf("You are possibly using a mesh/control file which was created before the removal.\n");
		printf("Please deprecate/delete any code functionality which allowed you to reach this point.\n");
		EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}

static void set_function_pointers_num_flux_bc_diffusion (struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_DIRICHLET: case BC_DIRICHLET_ALT1:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_diffusion_dirichlet;
		break;
	case BC_NEUMANN: case BC_NEUMANN_ALT1:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_diffusion_neumann;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}

static void set_function_pointers_num_flux_bc_euler (struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_RIEMANN:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_riemann;
		break;
	case BC_SLIPWALL:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_slipwall;
		break;
	case BC_SUPERSONIC_IN:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_supersonic_inflow;
		break;
	case BC_SUPERSONIC_OUT:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_supersonic_outflow;
		break;
	case BC_BACKPRESSURE:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_back_pressure;
		break;
	case BC_TOTAL_TP:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_euler_total_tp;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}

static void set_function_pointers_num_flux_bc_navier_stokes (struct Solver_Face_T* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	const int bc = face->bc % BC_STEP_SC;
	switch (bc) {
	case BC_RIEMANN:
	case BC_SUPERSONIC_IN:
	case BC_SUPERSONIC_OUT:
		set_function_pointers_num_flux_bc_euler(s_face);
		break;
	case BC_NOSLIP_ALL_ROTATING:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_navier_stokes_no_slip_all_rotating;
		break;
	case BC_NOSLIP_ADIABATIC:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_navier_stokes_no_slip_flux_adiabatic;
		break;
	case BC_NOSLIP_DIABATIC:
		s_face->constructor_Boundary_Value_fcl = constructor_Boundary_Value_T_navier_stokes_no_slip_flux_diabatic;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",face->bc);
		break;
	}
}

#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"

#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_boundary.h"
#include "undef_templates_test_case.h"
