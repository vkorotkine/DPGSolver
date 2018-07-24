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
#include "definitions_physics.h"
#include "definitions_tol.h"

#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "const_cast.h"
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
#include "solve_opg.h"

// Static function declarations ************************************************************************************* //

/** Flag for whether the solution should be reinitialized whenever \ref solve_for_solution is called.
 *
 *  \warning This should __only__ be enabled for testing purposes as it may significantly slow down convergence for
 *           non-linear problems if the initial solution is not close to the exact solution, which is true for all but
 *           the trivial cases.
 */
#define ALWAYS_SET_INITIAL false
/// \todo REMOVE THIS WHEN FINISHED TESTING! (Replace with restart file where applicable)

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
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solve_T.c"
#include "undef_templates_type.h"

void solve_for_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	if (ALWAYS_SET_INITIAL) {
		printf("*** Warning: Always resetting to initial solution. *** \n");
		set_initial_solution(sim);
#if 1
for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
	struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;
	set_to_value_Multiarray_d(s_vol->sol_coef,0.0);
	set_to_value_Multiarray_d(s_vol->test_s_coef,0.0);
}
#include "face_solver.h"
for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
	struct Solver_Face*const s_face = (struct Solver_Face*) curr;
	set_to_value_Multiarray_d(s_face->nf_coef,0.0);
}
#endif
	}

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

double compute_rlhs (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	double max_rhs = 0.0;

	zero_memory_volumes(sim->volumes);
	switch (sim->method) {
		case METHOD_DG:  max_rhs = compute_rlhs_dg(sim,ssi);    break;
		case METHOD_DPG: max_rhs = compute_rlhs_dpg(sim,ssi);   break;
		case METHOD_OPG: // fallthrough
		case METHOD_OPGC0:
			max_rhs = compute_rlhs_opg(sim,ssi);
			break;
		default:         EXIT_ERROR("Unsupported: %d\n",sim->method); break;
	}

	return max_rhs;
}

void copy_rhs (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	switch (sim->method) {
		case METHOD_DG: // fallthrough
		case METHOD_OPG: // fallthrough
		case METHOD_OPGC0:
			break; // Do nothing
		case METHOD_DPG: EXIT_ADD_SUPPORT; UNUSED(ssi); break;
		default:         EXIT_ERROR("Unsupported: %d\n",sim->method); break;
	}
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;

		destructor_conditional_Multiarray_d(s_vol->rhs_0);
		s_vol->rhs_0 = constructor_copy_Multiarray_d(s_vol->rhs); // keep
	}
}

void enforce_positivity_highorder (struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	if (!test_case_requires_positivity((struct Test_Case*) sim->test_case_rc->tc))
		return;

	struct Multiarray_d* s_coef_b = constructor_s_coef_bezier(s_vol,sim); // destructed

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

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	const struct Volume* vol         = (struct Volume*) s_vol;
	const struct Solver_Element* s_e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref;
	struct Multiarray_d* s_coef = s_vol->sol_coef;

	const struct Operator* ccBS0_vs_vs = get_Multiarray_Operator(s_e->ccBS0_vs_vs,(ptrdiff_t[]){0,0,p,p});
	mm_NN1C_Operator_Multiarray_d(
		ccBS0_vs_vs,(struct const_Multiarray_d*)s_coef_b,s_coef,op_format,s_coef_b->order,NULL,NULL);

	destructor_Multiarray_d(s_coef_b);
}

void destructor_Solver_Storage_Implicit (struct Solver_Storage_Implicit* ssi)
{
	if (!ssi->do_not_destruct_A)
		MatDestroy(&ssi->A);
	VecDestroy(&ssi->b);
	destructor_conditional_const_Vector_i(ssi->corr_l2_c0);
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
		dof += s_vol->sol_coef->extents[0];
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
	case 'l':
		dof += compute_dof_volumes_l_mult(sim);
		break;
	case 'k': // 'k'eep
		dof += compute_dof_faces(sim);
		dof += compute_dof_volumes_l_mult(sim);
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

double compute_max_rhs_from_ssi (const struct Solver_Storage_Implicit*const ssi)
{
	double max_rhs = 0.0;
	VecNorm(ssi->b,NORM_INFINITY,&max_rhs);
	return max_rhs;
}

void add_to_petsc_Mat (const struct Solver_Storage_Implicit*const ssi, const struct const_Matrix_d*const lhs)
{
	assert(lhs->layout == 'R');
	const ptrdiff_t ext_0 = lhs->ext_0,
	                ext_1 = lhs->ext_1;

	PetscInt idxm[ext_0],
	         idxn[ext_1];

	for (int i = 0; i < ext_0; ++i)
		idxm[i] = ssi->row+i;

	for (int i = 0; i < ext_1; ++i)
		idxn[i] = ssi->col+i;

	const PetscScalar*const vv = lhs->data;
	MatSetValues(ssi->A,(PetscInt)ext_0,idxm,(PetscInt)ext_1,idxn,vv,ADD_VALUES);
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
