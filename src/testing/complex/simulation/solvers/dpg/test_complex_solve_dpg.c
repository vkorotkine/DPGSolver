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

#include "test_complex_solve_dpg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_integration.h"

#include "test_complex_boundary.h"
#include "test_complex_boundary_advection.h"
#include "test_complex_compute_all_rhs_dpg.h"
#include "test_support_multiarray.h"

#include "face_solver_dpg_complex.h"
#include "volume_solver_dpg_complex.h"

#include "complex_multiarray.h"
#include "multiarray.h"

#include "boundary_advection.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Identical to \ref set_function_pointers_num_flux_dg for \ref Complex_DPG_Solver_Face\*s.
static void set_function_pointers_num_flux_dpg
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Compute the complex rhs terms for the dpg scheme.
static void compute_rhs_cmplx_step_dpg
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, ///< The \ref Complex_DPG_Solver_Volume.
	 struct Solver_Storage_Implicit* ssi,                 ///< \ref Solver_Storage_Implicit.
	 const struct Simulation* sim                         ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void perturb_solution_dpg (const struct Simulation* sim)
{
	assert(sim->test_case->has_2nd_order == false); // Add support.

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		perturb_Multiarray_d(s_vol->sol_coef,1e-5);
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face = (struct Solver_Face*) curr;
		perturb_Multiarray_d(s_face->nf_coef,1e-5);
	}
}

void set_initial_solution_complex_dpg (const struct Simulation* sim)
{
	assert(sim->test_case->has_2nd_order == false); // Add support.

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol                   = (struct Solver_Volume*) curr;
		struct Complex_DPG_Solver_Volume* c_dpg_s_vol = (struct Complex_DPG_Solver_Volume*) curr;

		set_Multiarray_c_Multiarray_d(c_dpg_s_vol->sol_coef,(struct const_Multiarray_d*)s_vol->sol_coef);
	}

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face                   = (struct Solver_Face*) curr;
		struct Complex_DPG_Solver_Face* c_dpg_s_face = (struct Complex_DPG_Solver_Face*) curr;

		set_Multiarray_c_Multiarray_d(c_dpg_s_face->nf_coef,(struct const_Multiarray_d*)s_face->nf_coef);
	}
}

void compute_lhs_cmplx_step_dpg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	set_function_pointers_num_flux_dpg(sim);

	for (struct Intrusive_Link* curr_c = sim->volumes->first; curr_c; curr_c = curr_c->next) {
		const struct Volume* vol                        = (struct Volume*) curr_c;
		struct Complex_DPG_Solver_Volume* c_dpg_s_vol_c = (struct Complex_DPG_Solver_Volume*) curr_c;
		struct Multiarray_c* sol_coef_c = c_dpg_s_vol_c->sol_coef;
		const ptrdiff_t n_col_l = compute_size(sol_coef_c->order,sol_coef_c->extents);
		for (int col_l = 0; col_l < n_col_l; ++col_l) {
			ssi->col = vol->index+col_l;

			sol_coef_c->data[col_l] += CX_STEP*I;
			compute_rhs_cmplx_step_dpg(c_dpg_s_vol_c,ssi,sim);
			sol_coef_c->data[col_l] -= CX_STEP*I;

//			set_col_lhs_cmplx_step_dg(col_l,(struct Solver_Volume*)curr_c,volumes_local,ssi);
EXIT_UNSUPPORTED;
		}
	}
	petsc_mat_vec_assemble(ssi);
// comment about directly storing imaginary part of the rhs in the ssi Vec. A different PETSc build affecting **all**
// containers would have to be used to support a complex Vector here.
EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void set_function_pointers_num_flux_dpg (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Solver_Face* s_face             = (struct Solver_Face*) curr;
		struct Complex_DPG_Solver_Face* c_dpg_s_face = (struct Complex_DPG_Solver_Face*) curr;

		if (s_face->constructor_Boundary_Value_fcl == constructor_Boundary_Value_s_fcl_interp)
			c_dpg_s_face->constructor_Boundary_Value_c_fcl = constructor_Boundary_Value_c_s_fcl_interp;
		else if (s_face->constructor_Boundary_Value_fcl == constructor_Boundary_Value_advection_inflow)
			c_dpg_s_face->constructor_Boundary_Value_c_fcl = constructor_Boundary_Value_c_advection_inflow;
		else if (s_face->constructor_Boundary_Value_fcl == constructor_Boundary_Value_advection_outflow)
			c_dpg_s_face->constructor_Boundary_Value_c_fcl = constructor_Boundary_Value_c_advection_outflow;
		else
			EXIT_UNSUPPORTED;
	}
}

static void compute_rhs_cmplx_step_dpg
	(const struct Complex_DPG_Solver_Volume* c_dpg_s_vol, struct Solver_Storage_Implicit* ssi,
	 const struct Simulation* sim)
{
	compute_all_rhs_dpg_c(c_dpg_s_vol,ssi,sim);
}
