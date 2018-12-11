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
#include "petscsys.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_integration_support_face.h"
#include "test_support.h"
#include "test_support_multiarray.h"
#include "test_support_solve.h"

#include "multiarray.h"

#include "face_solver.h"

#include "compute_face_rlhs.h"
#include "math_functions.h"
#include "numerical_flux.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the real \ref Numerical_Flux_T computed using the complex step method.
 *  \return See brief. */
static struct Numerical_Flux* constructor_Numerical_Flux_cmplx_step
	(const struct Solver_Face*const s_face_r,   ///< Real \ref Solver_Face_T.
	 const struct Solver_Face_c*const s_face_c, ///< Complex \ref Solver_Face_T.
	 struct Simulation* sim[2]                  ///< Real and complex \ref Simulation containers.
		);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the numerical fluxes (\ref test_integration_numerical_fluxes.cpp).
 *  \return 0 on success.
 *
 *  The linearizations are verified by comparing with the output when using the complex step method. Details of the
 *  complex step method can be found in Squire et al. \cite Squire1998 and Martins et al. \cite Martins2003.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 3,"Invalid number of input arguments");
	const char*const ctrl_name = argv[2];

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int p          = int_test_info->p_ref[0],
	          ml         = int_test_info->ml[0],
	          adapt_type = int_test_info->adapt_type;
	assert(adapt_type == ADAPT_0);

	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	const char type_rc[2] = { 'r', 'c', };
	struct Simulation* sim[2] = { NULL, NULL, };
	int ind_rc = -1;

	// real
	ind_rc = 0;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc],false); // destructed
	((struct Test_Case_c*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'i';
	struct Numerical_Flux_Input*const num_flux_i = constructor_Numerical_Flux_Input(sim[ind_rc]); // destructed
	perturb_solution(sim[ind_rc]);

	// complex
	ind_rc = 1;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc],false); // destructed
	((struct Test_Case_c*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'e';
	perturb_solution(sim[ind_rc]);

	bool pass_all = true;
	for (struct Intrusive_Link* curr_r = sim[0]->faces->first; curr_r; curr_r = curr_r->next) {
		struct Face* face_r          = (struct Face*) curr_r;
		struct Solver_Face* s_face_r = (struct Solver_Face*) curr_r;

		struct Face* face_c            = NULL;
		struct Solver_Face_c* s_face_c = NULL;
		for (struct Intrusive_Link* curr_c = sim[1]->faces->first; curr_c; curr_c = curr_c->next) {
			face_c = (struct Face*) curr_c;
			if (face_r->index == face_c->index) {
				s_face_c = (struct Solver_Face_c*) curr_c;
				break;
			}
		}

		constructor_Numerical_Flux_Input_data_with_gradients(num_flux_i,s_face_r,sim[0]); // destructed
		struct Numerical_Flux*const num_flux = constructor_Numerical_Flux(num_flux_i);    // destructed

		struct Numerical_Flux*const num_flux_c =
			constructor_Numerical_Flux_cmplx_step(s_face_r,s_face_c,sim); // destructed

		const struct Neigh_Info_NF*const ni   = num_flux->neigh_info;
		const struct Neigh_Info_NF*const ni_c = num_flux_c->neigh_info;

		const bool bdry      = face_r->boundary;
		const bool*const c_m = num_flux_i->flux_i->compute_member;

		assert(!c_m[3]); // Add support.
		assert(!c_m[4]); // Add support.
		assert(!c_m[5]); // Add support.
		bool pass        = false;
		const double tol = 3e1*EPS;
		/** \note The difference with no relative scaling is used below as some of the values were found to be
		 *        extremely small, resulting in round-off effects and requiring a large tolerance using the
		 *        standard difference check.
		 *
		 *  The necessity of doing this warrants further investigation/thought as removing the relative
		 *  contribution was not necessary when checking the full linearization. */
		const bool diffs[] =
			{ c_m[0]          ? diff_no_rel_const_Multiarray_d(num_flux->nnf,num_flux_c->nnf,tol) : false,
			  c_m[1]          ? diff_no_rel_const_Multiarray_d(ni[0].dnnf_ds,ni_c[0].dnnf_ds,tol) : false,
			  c_m[1] && !bdry ? diff_no_rel_const_Multiarray_d(ni[1].dnnf_ds,ni_c[1].dnnf_ds,tol) : false,
			  c_m[2]          ? diff_no_rel_const_Multiarray_d(ni[0].dnnf_dg,ni_c[0].dnnf_dg,tol) : false,
			  c_m[2] && !bdry ? diff_no_rel_const_Multiarray_d(ni[1].dnnf_dg,ni_c[1].dnnf_dg,tol) : false,
			};
		const int len = COUNT_OF(diffs);
		if (check_diff(len,diffs,&pass)) {
			if (diffs[0]) print_diff_no_rel_const_Multiarray_d(num_flux->nnf,num_flux_c->nnf,tol);
			if (diffs[1]) print_diff_no_rel_const_Multiarray_d(ni[0].dnnf_ds,ni_c[0].dnnf_ds,tol);
			if (diffs[2]) print_diff_no_rel_const_Multiarray_d(ni[1].dnnf_ds,ni_c[1].dnnf_ds,tol);
			if (diffs[3]) print_diff_no_rel_const_Multiarray_d(ni[0].dnnf_dg,ni_c[0].dnnf_dg,tol);
			if (diffs[4]) print_diff_no_rel_const_Multiarray_d(ni[1].dnnf_dg,ni_c[1].dnnf_dg,tol);
		}
		expect_condition(pass,"numerical_flux_linearization");
		if (!pass)
			pass_all = false;
#if 0
printf("info: %d %d %d %d %d %d\n",face_r->index,diffs[0],diffs[1],diffs[2],diffs[3],diffs[4]);
printf("info: %d %d \n",face_r->boundary,face_r->bc);
/* print_const_Multiarray_d(ni[0].dnnf_ds); */
/* print_const_Multiarray_d(ni_c[0].dnnf_ds); */
/* if (!bdry) print_diff_no_rel_const_Multiarray_d(ni[1].dnnf_ds,ni_c[1].dnnf_ds,tol); */
/* if (!pass) */
EXIT;
#endif

		destructor_Numerical_Flux_Input_data(num_flux_i);
		destructor_Numerical_Flux(num_flux);
		destructor_Numerical_Flux(num_flux_c);
	}
	assert_condition(pass_all);

	// complex
	ind_rc = 1;
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc],false);

	// real
	ind_rc = 0;
	destructor_Numerical_Flux_Input(num_flux_i);
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc],false);

	destructor_Integration_Test_Info(int_test_info);
	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Set the members of the input real \ref mutable_Numerical_Flux_T container to zero.
static void set_to_zero_Numerical_Flux
	(struct mutable_Numerical_Flux*const m_num_flux,    ///< Standard.
	 const struct Numerical_Flux_Input*const num_flux_i ///< Standard.
	);

static struct Numerical_Flux* constructor_Numerical_Flux_cmplx_step
	(const struct Solver_Face*const s_face_r, const struct Solver_Face_c*const s_face_c, struct Simulation* sim[2])
{
	UNUSED(s_face_c);
	struct Numerical_Flux_Input*const num_flux_i = constructor_Numerical_Flux_Input(sim[0]); // destructed
	constructor_Numerical_Flux_Input_data_with_gradients(num_flux_i,s_face_r,sim[0]);        // destructed

	struct mutable_Numerical_Flux*const m_num_flux =
	(struct mutable_Numerical_Flux*) constructor_Numerical_Flux(num_flux_i); // returned
	set_to_zero_Numerical_Flux(m_num_flux,num_flux_i);

	struct Numerical_Flux_Input_c*const num_flux_i_c = constructor_Numerical_Flux_Input_c(sim[1]); // destructed
	constructor_Numerical_Flux_Input_c_data_members(num_flux_i_c,num_flux_i,'b');                  // destructed

	const bool bdry      = ((struct Face*)s_face_r)->boundary;
	const bool*const c_m = num_flux_i->flux_i->compute_member;

	struct Multiarray_c*const s[] = { (struct Multiarray_c*) num_flux_i_c->bv_l.s,
	                                  !bdry ? (struct Multiarray_c*) num_flux_i_c->bv_r.s : NULL, };

	const ptrdiff_t n_n  = s[0]->extents[0];
	const int n_vr = get_set_n_var_eq(NULL)[0],
	          n_eq = get_set_n_var_eq(NULL)[1];

	for (int side = 0; side < 2; ++side) {
	for (int vr = 0; s[side] && vr < n_vr; ++vr) {
		double complex*const data_s = get_col_Multiarray_c(vr,s[side]);
		add_to_c(data_s,CX_STEP*I,n_n);
		if (bdry) { // Update right state if dependent on the left state.
			destructor_Boundary_Value_c_data(num_flux_i_c);
			constructor_Boundary_Value_c_data(num_flux_i_c,s_face_c,sim[1]); // destructed
		}
		struct Numerical_Flux_c* const num_flux_c = constructor_Numerical_Flux_c(num_flux_i_c); // destructed

		// nnf[NEQ]
		assert(c_m[0] == true);
		if (side == 0 && vr == 0) {
			const double complex*const data_c = get_col_const_Multiarray_c(0,num_flux_c->nnf);
			double*const data_r               = get_col_Multiarray_d(0,m_num_flux->nnf);
			for (int n = 0; n < n_n*n_eq; ++n)
				data_r[n] = creal(data_c[n]);
		}

		// dnnf_ds[NEQ,NVAR]
		if (c_m[1]) {
			struct m_Neigh_Info_NF*const ni_nf = &m_num_flux->neigh_info[side];

			for (int eq = 0; eq < n_eq; ++eq) {
				const ptrdiff_t ind_c = eq,
				                ind_r = eq+n_eq*(vr);

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,num_flux_c->nnf);
				double*const data_r               = get_col_Multiarray_d(ind_r,ni_nf->dnnf_ds);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}
		}
		destructor_Numerical_Flux_c(num_flux_c);
		add_to_c(data_s,-CX_STEP*I,n_n);
	}}

	struct Multiarray_c*const g[] = { (struct Multiarray_c*) num_flux_i_c->bv_l.g,
	                                  !bdry ? (struct Multiarray_c*) num_flux_i_c->bv_r.g : NULL, };

	if (c_m[2]) {
	for (int side = 0; side < 2; ++side) {
	for (int dx = 0; dx < DIM; ++dx) {
	for (int vr = 0; g[side] && vr < n_vr; ++vr) {
		double complex*const data_g = get_col_Multiarray_c(vr+n_vr*(dx),g[side]);
		add_to_c(data_g,CX_STEP*I,n_n);
		if (bdry) {
			destructor_Boundary_Value_c_data(num_flux_i_c);
			constructor_Boundary_Value_c_data(num_flux_i_c,s_face_c,sim[1]); // destructed
		}
		struct Numerical_Flux_c* const num_flux_c = constructor_Numerical_Flux_c(num_flux_i_c); // destructed

		struct m_Neigh_Info_NF*const ni_nf = &m_num_flux->neigh_info[side];

		// dnnf_dg[NEQ,NVAR,DIM]
		for (int eq = 0; eq < n_eq; ++eq) {
			const ptrdiff_t ind_c = eq,
			                ind_r = eq+n_eq*(vr+n_vr*(dx));

			const double complex*const data_c = get_col_const_Multiarray_c(ind_c,num_flux_c->nnf);
			double*const data_r               = get_col_Multiarray_d(ind_r,ni_nf->dnnf_dg);
			for (int n = 0; n < n_n; ++n)
				data_r[n] = cimag(data_c[n])/CX_STEP;
		}
		destructor_Numerical_Flux_c(num_flux_c);
		add_to_c(data_g,-CX_STEP*I,n_n);
	}}}}

	assert(!c_m[3]); // Add support.
	assert(!c_m[4]); // Add support.
	assert(!c_m[5]); // Add support.

	destructor_Numerical_Flux_Input_data(num_flux_i);
	destructor_Numerical_Flux_Input(num_flux_i);
	destructor_Numerical_Flux_Input_c_data_members(num_flux_i_c,'b');
	destructor_Numerical_Flux_Input_c(num_flux_i_c);

	return (struct Numerical_Flux*) m_num_flux;
}

// Level 1 ********************************************************************************************************** //

static void set_to_zero_Numerical_Flux
	(struct mutable_Numerical_Flux*const m_num_flux, const struct Numerical_Flux_Input*const num_flux_i)
{
	assert(sizeof(struct mutable_Numerical_Flux) == sizeof(struct Numerical_Flux));
	const bool*const c_m = num_flux_i->flux_i->compute_member;

	if (c_m[0])
		set_to_value_Multiarray_d(m_num_flux->nnf,0.0);
	if (c_m[1]) {
		set_to_value_Multiarray_d(m_num_flux->neigh_info[0].dnnf_ds,0.0);
		if (m_num_flux->neigh_info[1].dnnf_ds)
			set_to_value_Multiarray_d(m_num_flux->neigh_info[1].dnnf_ds,0.0);
	}
	if (c_m[2]) {
		set_to_value_Multiarray_d(m_num_flux->neigh_info[0].dnnf_dg,0.0);
		if (m_num_flux->neigh_info[1].dnnf_ds)
			set_to_value_Multiarray_d(m_num_flux->neigh_info[1].dnnf_dg,0.0);
	}
	if (c_m[3])
		EXIT_ADD_SUPPORT;
	if (c_m[4])
		EXIT_ADD_SUPPORT;
	if (c_m[5])
		EXIT_ADD_SUPPORT;
}
