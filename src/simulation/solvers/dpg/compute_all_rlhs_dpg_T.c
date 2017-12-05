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
 *  \todo update includes.
 *  \todo Attempt to remove all preprocessor conditionals.
 */

// Static function declarations ************************************************************************************* //

/// \brief Increment and add dof for the rhs and lhs with the face contributions from 1st order equations.
static void increment_and_add_dof_rlhs_f_1_T
	(struct Vector_T* rhs,                      ///< Holds the values of the rhs.
	 struct Matrix_T** lhs_ptr,                 ///< Pointer to the matrix holding the values of the lhs.
	 const struct DPG_Solver_Volume* dpg_s_vol, ///< \ref DPG_Solver_Volume.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Increment the rhs and lhs entries corresponding to internal faces.
static void increment_rlhs_internal_face_T
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< The current \ref DPG_Solver_Volume.
	 const struct DPG_Solver_Face* dpg_s_face,  ///< The current \ref DPG_Solver_Face.
	 struct Matrix_T* lhs,                      ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_T* rhs,                      ///< The rhs matrix contribution for the current volume/faces.
	 int* ind_dof,                              ///< The index of the current dof under consideration.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

/// \brief Increment the rhs and lhs entries corresponding to boundary faces.
static void increment_rlhs_boundary_face_T
	(const struct DPG_Solver_Volume* dpg_s_vol, ///< The current \ref DPG_Solver_Volume.
	 const struct DPG_Solver_Face* dpg_s_face,  ///< The current \ref DPG_Solver_Face.
	 struct Matrix_T* lhs,                      ///< The lhs matrix contribution for the current volume/faces.
	 struct Matrix_T* rhs,                      ///< The rhs matrix contribution for the current volume/faces.
	 const struct Simulation* sim               ///< \ref Simulation.
	);

static void increment_and_add_dof_rlhs_f_1_T
	(struct Vector_T* rhs, struct Matrix_T** lhs_ptr, const struct DPG_Solver_Volume* dpg_s_vol,
	 const struct Simulation* sim)
{
	struct Matrix_T* lhs = *lhs_ptr;

	const struct Volume* vol          = (struct Volume*) dpg_s_vol;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) dpg_s_vol;

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	const ptrdiff_t n_dof_s  = (lhs->ext_1)/n_vr,
	                n_dof_nf = compute_n_dof_nf(s_vol);
#if TYPE_RC == TYPE_REAL
	struct Matrix_d* lhs_add = constructor_empty_Matrix_d('R',lhs->ext_0,(n_dof_s+n_dof_nf)*n_vr); // moved
	set_to_value_Matrix_d(lhs_add,0.0);
#elif TYPE_RC == TYPE_COMPLEX
	struct Matrix_c* lhs_add = constructor_empty_Matrix_c('R',lhs->ext_0,(n_dof_s+n_dof_nf)*n_vr); // moved
#endif
	set_block_Matrix_T(lhs_add,(struct const_Matrix_T*)lhs,0,0,'i');
	destructor_Matrix_T(lhs);
	lhs = lhs_add;

	struct Matrix_T rhs_M =
		{ .layout = 'C', .ext_0 = (rhs->ext_0)/n_eq, .ext_1 = n_eq, .owns_data = false, .data = rhs->data, };

	int ind_dof = (int)(n_vr*n_dof_s);
	for (int i = 0; i < NFMAX;    ++i) {
	for (int j = 0; j < NSUBFMAX; ++j) {
		const struct Face* face = vol->faces[i][j];
		if (!face)
			continue;

		const struct DPG_Solver_Face* dpg_s_face = (struct DPG_Solver_Face*) face;
		if (!face->boundary)
			increment_rlhs_internal_face_T(dpg_s_vol,dpg_s_face,lhs,&rhs_M,&ind_dof,sim);
		else
			increment_rlhs_boundary_face_T(dpg_s_vol,dpg_s_face,lhs,&rhs_M,sim);
	}}
//print_Matrix_T(lhs);
	*lhs_ptr = lhs;
}

// Level 1 ********************************************************************************************************** //

/// \brief Scale required \ref Numerical_Flux terms by the face Jacobian.
static void scale_by_Jacobian_T
	(const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face         ///< The current \ref Solver_Face.
	);

/// \brief Increment the rhs terms with the contribution of the boundary face.
static void increment_rhs_boundary_face_T
	(struct Matrix_T* rhs,                    ///< Holds the rhs terms.
	 const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face,        ///< The current \ref Solver_Face.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

/// \brief Increment the lhs terms with the contribution of the boundary face.
static void increment_lhs_boundary_face_T
	(struct Matrix_T* lhs,                    ///< Holds the lhs terms.
	 const struct Numerical_Flux_T* num_flux, ///< \ref Numerical_Flux.
	 const struct Solver_Face* s_face,        ///< The current \ref Solver_Face.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

static void increment_rlhs_internal_face_T
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct DPG_Solver_Face* dpg_s_face, struct Matrix_T* lhs,
	 struct Matrix_T* rhs, int* ind_dof, const struct Simulation* sim)
{
	/// As the rhs is **always** linear wrt the trace unknowns, the rhs and lhs are computed together.
	const struct const_Matrix_R* lhs_l = constructor_lhs_l_internal_face_dpg(dpg_s_vol,dpg_s_face); // destructed
//print_const_Matrix_d(tw0_vt_fc_op->op_std);
#if TYPE_RC == TYPE_COMPLEX
	const struct Complex_DPG_Solver_Face* s_face = (struct Complex_DPG_Solver_Face*) dpg_s_face;
#elif TYPE_RC == TYPE_REAL
	const struct Solver_Face* s_face  = (struct Solver_Face*) dpg_s_face;
#endif
	struct Matrix_T nf_coef = interpret_Multiarray_as_Matrix_T(s_face->nf_coef);
	mm_RTT('N','N',1.0,1.0,lhs_l,(struct const_Matrix_T*)&nf_coef,rhs);

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;

	const ptrdiff_t n_dof_test = (lhs->ext_0)/n_eq,
	                n_dof_nf   = nf_coef.ext_0;
	for (int vr = 0; vr < n_vr; ++vr) {
		set_block_Matrix_T_R(lhs,lhs_l,vr*n_dof_test,*ind_dof,'i');
		*ind_dof += (int)n_dof_nf;
	}
	destructor_const_Matrix_d(lhs_l);
}

static void increment_rlhs_boundary_face_T
	(const struct DPG_Solver_Volume* dpg_s_vol, const struct DPG_Solver_Face* dpg_s_face, struct Matrix_T* lhs,
	 struct Matrix_T* rhs, const struct Simulation* sim)
{
	UNUSED(dpg_s_vol);
	struct Numerical_Flux_Input_T* num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed

	const struct Solver_Face* s_face = (struct Solver_Face*) dpg_s_face;
#if TYPE_RC == TYPE_COMPLEX
	const struct Complex_DPG_Solver_Face* c_dpg_s_face = (struct Complex_DPG_Solver_Face*) dpg_s_face;
	constructor_Numerical_Flux_Input_c_data_dpg(num_flux_i,c_dpg_s_face,sim); // destructed
#elif TYPE_RC == TYPE_REAL
	constructor_Numerical_Flux_Input_data(num_flux_i,s_face,sim); // destructed
#endif
	struct Numerical_Flux_T* num_flux = constructor_Numerical_Flux_T(num_flux_i); // destructed
	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	scale_by_Jacobian_T(num_flux,s_face);

	increment_rhs_boundary_face_T(rhs,num_flux,s_face,sim);
	increment_lhs_boundary_face_T(lhs,num_flux,s_face,sim);

	destructor_Numerical_Flux_T(num_flux);
}

// Level 2 ********************************************************************************************************** //

static void scale_by_Jacobian_T (const struct Numerical_Flux_T* num_flux, const struct Solver_Face* s_face)
{
	const struct Face* face = (struct Face*) s_face;

	assert(face->boundary);
	assert(num_flux->neigh_info[0].dnnf_ds != NULL || num_flux->neigh_info[0].dnnf_dg != NULL);

	const struct const_Vector_R jacobian_det_fc = interpret_const_Multiarray_as_Vector_R(s_face->jacobian_det_fc);
	scale_Multiarray_T_by_Vector_R('L',1.0,(struct Multiarray_T*)num_flux->nnf,&jacobian_det_fc,false);

	if (num_flux->neigh_info[0].dnnf_ds)
		scale_Multiarray_T_by_Vector_R(
			'L',1.0,(struct Multiarray_T*)num_flux->neigh_info[0].dnnf_ds,&jacobian_det_fc,false);
	if (num_flux->neigh_info[0].dnnf_dg)
		EXIT_ADD_SUPPORT;
}

static void increment_rhs_boundary_face_T
	(struct Matrix_T* rhs, const struct Numerical_Flux_T* num_flux, const struct Solver_Face* s_face,
	 const struct Simulation* sim)
{
	ptrdiff_t extents[2] = { rhs->ext_0, rhs->ext_1, };
	struct Multiarray_T rhs_Ma =
		{ .layout = 'C', .order = 2, .extents = extents, .owns_data = false, .data = rhs->data, };

	const struct Operator* tw0_vt_fc = get_operator__tw0_vt_fc(0,s_face);

	UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';
	mm_NNC_Operator_Multiarray_T(-1.0,1.0,tw0_vt_fc,num_flux->nnf,&rhs_Ma,op_format,2,NULL,NULL);
}

static void increment_lhs_boundary_face_T
	(struct Matrix_T* lhs, const struct Numerical_Flux_T* num_flux, const struct Solver_Face* s_face,
	 const struct Simulation* sim)
{
	UNUSED(sim);
	assert(((struct Face*)s_face)->boundary);

	struct Matrix_T* lhs_ll = constructor_lhs_f_1_T((int[]){0,0},num_flux,s_face); // destructed

	set_block_Matrix_T(lhs,(struct const_Matrix_T*)lhs_ll,0,0,'a');
#if 0
printf("lhs\n");
print_Matrix_T(lhs_ll);
print_Matrix_T(lhs);
#endif
	destructor_Matrix_T(lhs_ll);
}
