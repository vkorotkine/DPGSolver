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

#include "solve.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "geometry.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_explicit.h"
#include "solve_implicit.h"
#include "solution.h"
#include "solution_euler.h"
#include "test_case.h"

#include "solve_dg.h"
#include "solve_dpg.h"

// Static function declarations ************************************************************************************* //

/// \brief Set the memory of the rhs and lhs (if applicable) terms to zero for the volumes.
static void zero_memory_volumes
	(struct Intrusive_List* volumes ///< The list of volumes for which to set the memory.
	);

/// \brief Set the memory of the rhs and lhs (if applicable) terms relating to flux imbalances to zero for the volumes.
static void zero_memory_volumes_flux_imbalances
	(struct Intrusive_List* volumes ///< The list of volumes for which to set the memory.
	);

/// \brief Correct the coefficients such that `coef = t*coef + (1-t)*coef_avg`.
static void correct_coef
	(struct Multiarray_d*const coef, ///< The multiarray of coefficients for each of the variables.
	 const double p_ho               ///< The percentage of the 'h'igh-'o'rder contribution to retain.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_T.c"

void solve_for_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	set_up_solver_geometry(sim);
	set_initial_solution(sim);

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	switch (test_case->solver_proc) {
	case SOLVER_E:
		solve_explicit(sim);
		break;
	case SOLVER_I:
		solve_implicit(sim);
		break;
	case SOLVER_EI:
		solve_explicit(sim);
		solve_implicit(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
}

double compute_rhs (const struct Simulation* sim)
{
/// \todo Add assertions relevant to rhs.
	double max_rhs = 0.0;

	zero_memory_volumes(sim->volumes);
	switch (sim->method) {
	case METHOD_DG:
		max_rhs = compute_rhs_dg(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	return max_rhs;
}

double compute_rlhs (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	double max_rhs = 0.0;

	zero_memory_volumes(sim->volumes);
	switch (sim->method) {
		case METHOD_DG:  max_rhs = compute_rlhs_dg(sim,s_store_i);    break;
		case METHOD_DPG: max_rhs = compute_rlhs_dpg(sim,s_store_i);   break;
		default:         EXIT_ERROR("Unsupported: %d\n",sim->method); break;
	}

	return max_rhs;
}

void enforce_positivity_highorder (struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	if (!test_case_requires_positivity((struct Test_Case*) sim->test_case_rc->tc))
		return;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	const struct Volume* vol         = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref;
	const struct Operator* ccSB0_vs_vs = get_Multiarray_Operator(s_e->ccSB0_vs_vs,(ptrdiff_t[]){0,0,p,p});

	struct Multiarray_d* s_coef = s_vol->sol_coef;

	struct Multiarray_d* s_coef_b =
		constructor_mm_NN1_Operator_Multiarray_d(ccSB0_vs_vs,s_coef,'C',op_format,s_coef->order,NULL); // dest.

	convert_variables(s_coef_b,'c','p');

	// Note: As the Bezier basis forms a partition of unity, the average value is given by the average of the
	//       coefficients.
	const ptrdiff_t n_n  = s_coef_b->extents[0],
	                n_vr = s_coef_b->extents[1];

	for (int vr = 0; vr < n_vr; vr += (int)n_vr-1) {
		double* vr_data = &s_coef_b->data[vr*n_n];

		double vr_avg = average_d(vr_data,n_n);
		if (vr_avg < EPS_PHYS) {
			print_Multiarray_d(s_coef_b);
			EXIT_ERROR("Average %s approaching 0.0.\n",(vr == 0 ? "density" : "pressure"));
		}

		// Correct terms if necessary
		const double vr_min = minimum_d(vr_data,n_n);
		if (vr_min < EPS_PHYS) {
			const double p_ho = GSL_MAX(0.0,GSL_MIN(1.0,(vr_avg-EPS_PHYS)/(vr_avg-vr_min)));
			correct_coef(s_coef_b,p_ho);
		}
	}
	convert_variables(s_coef_b,'p','c');

	const struct Operator* ccBS0_vs_vs = get_Multiarray_Operator(s_e->ccBS0_vs_vs,(ptrdiff_t[]){0,0,p,p});
	mm_NN1C_Operator_Multiarray_d(
		ccBS0_vs_vs,(struct const_Multiarray_d*)s_coef_b,s_coef,op_format,s_coef_b->order,NULL,NULL);

	destructor_Multiarray_d(s_coef_b);
}

void destructor_Solver_Storage_Implicit (struct Solver_Storage_Implicit* ssi)
{
	MatDestroy(&ssi->A);
	VecDestroy(&ssi->b);

	free(ssi);
}

void increment_nnz (struct Vector_i* nnz, const ptrdiff_t ind_dof, const ptrdiff_t n_row, const ptrdiff_t n_col)
{
	assert(ind_dof >= 0);

	const ptrdiff_t i_max = ind_dof+n_row;
	assert(i_max <= nnz->ext_0);

	for (ptrdiff_t i = ind_dof; i < i_max; ++i)
		nnz->data[i] += (int)n_col;
}

void petsc_mat_vec_assemble (struct Solver_Storage_Implicit* ssi)
{
	MatAssemblyBegin(ssi->A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(ssi->A,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(ssi->b);
	VecAssemblyEnd(ssi->b);
}

ptrdiff_t compute_dof_sol_1st (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		dof += compute_size(s_vol->sol_coef->order,s_vol->sol_coef->extents);
	}
	return dof;
}

ptrdiff_t compute_dof_schur (const char dof_type, const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	switch (dof_type) {
	case 'f':
		dof += compute_dof_faces(sim);
		break;
	case 'v':
		dof += compute_dof_volumes(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",dof_type);
		break;
	}
	return dof;
}

void compute_flux_imbalances (struct Simulation*const sim)
{
	zero_memory_volumes_flux_imbalances(sim->volumes);

	switch (sim->method) {
	case METHOD_DG:
		compute_flux_imbalances_dg(sim);
		break;
	case METHOD_DPG:
		compute_flux_imbalances_dpg(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void zero_memory_volumes (struct Intrusive_List* volumes)
{
	zero_memory_volumes_flux_imbalances(volumes);
}

static void zero_memory_volumes_flux_imbalances (struct Intrusive_List* volumes)
{
	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
		set_to_value_Vector_d(s_vol->flux_imbalance,0.0);
	}
}

static void correct_coef (struct Multiarray_d*const coef, const double p_ho)
{
	assert(coef->order == 2);

	const ptrdiff_t n_n  = coef->extents[0],
	                n_vr = coef->extents[1];

	for (int vr = 0; vr < n_vr; ++vr) {
		double*const data = &coef->data[vr*n_n];
		double avg = average_d(data,n_n);
		if ((vr == 0 || vr == n_vr-1) && avg < EPS_PHYS) {
			assert(p_ho < EPS);
			avg = EPS_PHYS;
		}

		for (int n = 0; n < n_n; ++n)
			data[n] = p_ho*data[n] + (1-p_ho)*avg;
	}
}
