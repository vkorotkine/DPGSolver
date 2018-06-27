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

#include "solve_opg.h"

#include <assert.h>

#include "macros.h"
#include "definitions_test_case.h"

#include "face_solver.h"
#include "volume_solver_opg.h"

#include "multiarray.h"
#include "vector.h"

#include "compute_face_rlhs_opg.h"
#include "compute_rlhs.h"
#include "compute_volume_rlhs_opg.h"
#include "const_cast.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Fill \ref Solver_Storage_Implicit::b with the negated rhs values.
 *
 *  See comments in \ref solve_implicit.h for why the negated values are used here.
 */
static void fill_petsc_Vec_b_opg
	(const struct Simulation*const sim,       ///< Standard.
	 struct Solver_Storage_Implicit*const ssi ///< Standard.
	);

/** \brief Convert \ref Solver_Storage_Implicit::A and Solver_Storage_Implicit::b to correspond to the globally
 *         C0 continuous dof. */
static void convert_petsc_to_c0_opg
	(struct Solver_Storage_Implicit*const ssi, ///< Standard.
	 const struct Simulation*const sim         ///< Standard.
	 );

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_opg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solve_opg_T.c"
#include "undef_templates_type.h"

double compute_rlhs_opg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	initialize_zero_memory_volumes(sim->volumes);
	compute_volume_rlhs_opg(sim,ssi,sim->volumes);
	switch (get_set_method(NULL)) {
	case METHOD_OPG:
		compute_face_rlhs_opg(sim,ssi,sim->faces);
		break;
	case METHOD_OPGC0:
		compute_face_rlhs_opg_boundary_d(sim,ssi,sim->faces);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL));
		break;
	}

	compute_source_rhs_dg_like(sim);
	fill_petsc_Vec_b_opg(sim,ssi);

	struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;
	switch (test_case->lhs_terms) {
	case LHS_FULL_NEWTON:
		break; // Do nothing
	case LHS_CFL_RAMPING:
		printf("*** Warning: CFL ramping not performed for the OPG method. ***\n");
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->lhs_terms);
		break;
	}
	const double max_rhs = compute_max_rhs_dg_like(sim);

	switch (get_set_method(NULL)) {
	case METHOD_OPG:
		break; // do nothing.
	case METHOD_OPGC0:
		convert_petsc_to_c0_opg(ssi,sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL));
		break;
	}

	return max_rhs;
}

void set_petsc_Mat_row_col_opg
	(struct Solver_Storage_Implicit*const ssi, const struct OPG_Solver_Volume*const v_l, const int eq,
	 const struct OPG_Solver_Volume*const v_r, const int vr)
{
	const struct Solver_Volume*const sv_l = (struct Solver_Volume*) v_l,
	                          *const sv_r = (struct Solver_Volume*) v_r;
	ssi->row = (int)(sv_l->ind_dof_test+sv_l->test_s_coef->extents[0]*eq);
	ssi->col = (int)(sv_r->ind_dof_test+sv_r->test_s_coef->extents[0]*vr);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the container of the redundant xyz coordinates corresponding to the L2 degrees of freedom
 *         with an additional index corresponding to each individual node.
 *  \return See brief. */
static const struct const_Multiarray_Vector_d* constructor_xyz_p_ind_L2
	(const struct Simulation*const sim ///< Standard.
	 );

static void fill_petsc_Vec_b_opg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const int ind_dof        = (int)((struct Solver_Volume*)curr)->ind_dof_test;
		const struct Multiarray_d*const rhs = ((struct Solver_Volume*)curr)->rhs;

		const int ni = (int)compute_size(rhs->order,rhs->extents);

		PetscInt    ix[ni];
		PetscScalar y[ni];

		for (int i = 0; i < ni; ++i) {
			ix[i] = ind_dof+i;
			y[i]  = -(rhs->data[i]);
		}

		VecSetValues(ssi->b,ni,ix,y,INSERT_VALUES);
	}
}

static void convert_petsc_to_c0_opg (struct Solver_Storage_Implicit*const ssi, const struct Simulation*const sim)
{
	assert((strcmp(sim->basis_geom,"lagrange") == 0));
	const struct const_Multiarray_Vector_d*const xyz_vt_red = constructor_xyz_p_ind_L2(sim);
	print_const_Multiarray_Vector_d(xyz_vt_red);

	// 1. Constructor (Vector_i) inds unique (ind of last entry is total number of unique entries).
	// 2. Constructor correspondence L2 to C0 (Vector_i: ext_0 = dof L2); store in ssi.
	// 3. Constructor for petsc Mat to transform L2 to C0.

	EXIT_ADD_SUPPORT; UNUSED(ssi);
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_Vector_d* constructor_xyz_p_ind_L2 (const struct Simulation*const sim)
{
	struct Vector_d*const xyz_vt_data = constructor_empty_Vector_d(0); // destructed.
	ptrdiff_t n_dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume*const s_vol         = (struct Solver_Volume*) curr;
		const struct OPG_Solver_Volume*const opg_s_vol = (struct OPG_Solver_Volume*) curr;
		const struct Operator*const cv0_vg_vt = get_operator__cv0_vg_vt_d(opg_s_vol);

		const char op_format = get_set_op_format(0);
		const struct const_Multiarray_d*const xyz_vt =
			constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vg_vt,s_vol->geom_coef,'R',op_format,2,NULL); // d.

		const ptrdiff_t n_vt = xyz_vt->extents[0];
		struct Vector_d*const xyz_vt_p_ind_V = constructor_empty_Vector_d(n_vt*(DIM+1)); // destructed
		for (int i = 0, ind = 0; i < n_vt; ++i) {
			for (int d = 0; d < DIM; ++d)
				xyz_vt_p_ind_V->data[ind++] = get_row_const_Multiarray_d(i,xyz_vt)[d];
			xyz_vt_p_ind_V->data[ind++] = (double) n_dof+i;
		}

		push_back_Vector_Vector_d(xyz_vt_data,(struct const_Vector_d*)xyz_vt_p_ind_V);
		n_dof += n_vt;

		destructor_const_Multiarray_d(xyz_vt);
		destructor_Vector_d(xyz_vt_p_ind_V);
	}
	struct Vector_i*const ext_V = constructor_empty_Vector_i(n_dof);
	set_to_value_Vector_i(ext_V,DIM+1);

	struct Multiarray_Vector_d*const xyz_p_ind = constructor_empty_Multiarray_Vector_d(true,1,(ptrdiff_t[]){n_dof}); // r.
	set_Multiarray_Vector_d_d(xyz_p_ind,xyz_vt_data->data,ext_V->data);
	destructor_Vector_d(xyz_vt_data);

	sort_Multiarray_Vector_tol_d(xyz_p_ind,false,1e-10);

	return (struct const_Multiarray_Vector_d*) xyz_p_ind;
}
