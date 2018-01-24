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
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_test_case.h"

#include "def_templates_compute_grad_coef_dg.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_face_solver.h"
#include "def_templates_face_solver_dg.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

#include "def_templates_boundary.h"
#include "def_templates_compute_face_rlhs.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solve.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Compute \ref DG_Solver_Volume_T::grad_coef_v for the list of volumes.
static void compute_grad_coef_volumes
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Compute \ref DG_Solver_Face_T::grad_coef_f for the list of faces.
static void compute_grad_coef_faces
	(const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_grad_coef_dg_T (const struct Simulation*const sim)
{
	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	if (!test_case->has_2nd_order)
		return;

	compute_grad_coef_volumes(sim);
	compute_grad_coef_faces(sim);

	EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the appropriate sub-range of the \ref DG_Solver_Element::cv1_vs_vc operators.
 *  \return See brief. */
static struct Multiarray_Operator get_operator__cv1_vs_vc
	(const struct DG_Solver_Volume_T*const dg_s_vol ///< \ref DG_Solver_Volume_T.
	);

/** \brief Constructor for the partial physical gradient operator from the reference operator and metric terms.
 *  \return See brief.
 *
 *  A partial operator is returned in the sense that the inverse Jacobian determinant contribution is omitted.
 */
static const struct const_Matrix_R* constructor_grad_xyz_p
	(const int dir,                                   ///< The direction. Options: 0 (x), 1 (y), 2 (z).
	 const struct Multiarray_Operator*const grad_rst, ///< The reference gradient operators.
	 const struct const_Multiarray_R*const metrics    ///< The metric terms.
	);

/** \brief Constructor for the \ref DG_Solver_Face_T d_g_coef_f__d_s_coef term associated with the current side and
 *         solution indices for 'i'nternal faces.
 *
 *  The `side_index` input determines the \ref DG_Solver_Face::neigh_info, while the `sol_index` input
 *  determines the \ref Neigh_Info_DG::d_g_coef_f__d_s_coef.
 */
static void constructor_d_g_coef_f__d_s_coef_i
	(const int side_index,                        ///< Defined for \ref get_sol_scale.
	 const int sol_index,                         ///< Defined for \ref get_sol_scale.
	 const int ind_num_flux_2nd,                  ///< Defined for \ref get_sol_scale.
	 const struct const_Matrix_R*const jdet_n_fc, ///< Jacobian determinant dotted with normals at fc nodes.
	 struct DG_Solver_Face_T*const dg_s_face,     ///< \ref DG_Solver_Face_T.
	 const bool collocated                        ///< \ref Simulation::collocated.
	);

/// \brief Constructor for the \ref DG_Solver_Face_T d_g_coef_f__d_s_coef term for 'b'oundary faces.
static void constructor_d_g_coef_f__d_s_coef_b
	(const int ind_num_flux_2nd,                  ///< Defined for \ref constructor_d_g_coef_f__d_s_coef_i.
	 const struct const_Matrix_R*const jdet_n_fc, ///< Defined for \ref constructor_d_g_coef_f__d_s_coef_i.
	 struct DG_Solver_Face_T*const dg_s_face,     ///< Defined for \ref constructor_d_g_coef_f__d_s_coef_i.
	 const struct Simulation*const sim            ///< \ref Simulation.
	);

static void compute_grad_coef_volumes (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) curr;
		struct DG_Solver_Volume_T*const dg_s_vol = (struct DG_Solver_Volume_T*) curr;

		for (int d = 0; d < DIM; ++d) {
			const struct Multiarray_Operator cv1_vs_vc = get_operator__cv1_vs_vc(dg_s_vol);
			const struct const_Matrix_R*const grad_xyz =
				constructor_grad_xyz_p(d,&cv1_vs_vc,s_vol->metrics_vc); // destructed/keep

			const struct const_Matrix_R* d_g_coef_v__d_s_coef = NULL;
			if (!sim->collocated) {
				const struct Operator*const tw0_vt_vc = get_operator__tw0_vt_vc_T(s_vol);
				const struct const_Matrix_R*const ibp2_v =
					constructor_mm_const_Matrix_R('N','N',1.0,tw0_vt_vc->op_std,grad_xyz,'R'); // destructed
				d_g_coef_v__d_s_coef =
					constructor_mm_const_Matrix_R('N','N',1.0,dg_s_vol->m_inv,ibp2_v,'R'); // keep
				destructor_const_Matrix_R(ibp2_v);
			} else {
				// Multiplication by the cubature weights is **not** performed here as multiplication by the
				// inverse cubature weights is **not** performed when multiplying by the inverse mass matrix.
				const struct const_Vector_R jacobian_det_vc =
					interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);
				scale_Matrix_R_by_Vector_R('L',1.0,(struct Matrix_R*)grad_xyz,&jacobian_det_vc,true);
				d_g_coef_v__d_s_coef = grad_xyz;
			}

			struct Multiarray_T grad_coef_v =
				interpret_Multiarray_as_slice_T(dg_s_vol->grad_coef_v,2,(ptrdiff_t[]){d});
			mm_NNC_Multiarray_T(1.0,0.0,d_g_coef_v__d_s_coef,
			                    (struct const_Multiarray_T*)s_vol->sol_coef,&grad_coef_v);
			copy_into_Multiarray_T(s_vol->grad_coef,(struct const_Multiarray_T*)dg_s_vol->grad_coef_v);

			const struct Test_Case_T*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
			if (test_case->solver_method_curr == 'i') {
				destructor_const_Matrix_T(dg_s_vol->d_g_coef_v__d_s_coef[d]);
				dg_s_vol->d_g_coef_v__d_s_coef[d] = d_g_coef_v__d_s_coef;
				if (!sim->collocated)
					destructor_const_Matrix_R(grad_xyz);
			} else {
				assert(test_case->solver_method_curr == 'e');
				destructor_const_Matrix_R(grad_xyz);
			}
		}
	}
}

static void compute_grad_coef_faces (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face*const face                  = (struct Face*) curr;
		struct Solver_Face_T*const s_face       = (struct Solver_Face_T*) curr;
		struct DG_Solver_Face_T*const dg_s_face = (struct DG_Solver_Face_T*) curr;

		const struct const_Vector_R jdet_fc    = interpret_const_Multiarray_as_Vector_R(s_face->jacobian_det_fc);
		const struct const_Matrix_R normals_fc = interpret_const_Multiarray_as_Matrix_R(s_face->normals_fc);

		struct Matrix_R*const jdet_n_fc = constructor_mm_diag_Matrix_R(1.0,&normals_fc,&jdet_fc,'L',false); // dest.
		if (jdet_n_fc->layout != 'C')
			transpose_Matrix_R(jdet_n_fc,true);

		const struct Test_Case_T*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
		const int ind_num_flux = test_case->ind_num_flux[1];
		if (!face->boundary) {
			const bool coll = sim->collocated;
			constructor_d_g_coef_f__d_s_coef_i(0,0,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);
			constructor_d_g_coef_f__d_s_coef_i(0,1,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);

			// Permute rows and negate such that the normals are interpreted as having the opposite volume be
			// dominant.
			transpose_Matrix_R(jdet_n_fc,true);
			assert(jdet_n_fc->layout == 'R');
			permute_Matrix_T_V(jdet_n_fc,get_operator__nc_fc_T(0,s_face));
			transpose_Matrix_R(jdet_n_fc,true);

			scale_Matrix_R(jdet_n_fc,-1.0);

			constructor_d_g_coef_f__d_s_coef_i(1,0,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);
			constructor_d_g_coef_f__d_s_coef_i(1,1,ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,coll);

EXIT_ADD_SUPPORT; // Compute the grad_coef_f terms (note linear in left/right solutions).
		} else {
			constructor_d_g_coef_f__d_s_coef_b(ind_num_flux,(struct const_Matrix_R*)jdet_n_fc,dg_s_face,sim);
EXIT_ADD_SUPPORT; // Compute the grad_coef_f term (likely as part of the function above).
		}

		destructor_Matrix_R(jdet_n_fc);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for a matrix holding the values of \ref DG_Solver_Volume_T::m_inv multiplied by the appropriate
 *         \ref Solver_Element::tw0_vt_fc operator.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_m_inv_tw0_vt_fc
	(const int side_index,                    ///< The side index of the volume associated with the face.
	 struct DG_Solver_Face_T*const dg_s_face, ///< \ref DG_Solver_Face_T.
	 const bool collocated                    ///< \ref Simulation::collocated.
	);

/** \brief Return a statically allocated array of \ref const_Vector_R\*s with each holding the data in one column of the
 *         input matrix.
 *  \return See brief. */
static struct const_Vector_R* get_jn_fc_V
	(const struct const_Matrix_R*const jdet_n_fc ///< The input matrix.
	);

/** \brief Compute the scaling associated with the \ref Solver_Volume_T::sol_coef contribution to the current face weak
 *         gradient term in the integrated by parts twice form. */
static Real get_sol_scale
	(const int side_index,      ///< The side index of the volume associated with the face.
	 const int sol_index,       ///< The side index of the solution.
	 const int ind_num_flux_2nd ///< The second component of \ref Test_Case_T::ind_num_flux.
	);

static struct Multiarray_Operator get_operator__cv1_vs_vc (const struct DG_Solver_Volume_T*const dg_s_vol)
{
	const struct Volume*const vol               = (struct Volume*) dg_s_vol;
	const struct Solver_Volume_T*const s_vol    = (struct Solver_Volume_T*) dg_s_vol;
	const struct DG_Solver_Element*const dg_s_e = (struct DG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;

	return set_MO_from_MO(dg_s_e->cv1_vs_vc[curved],1,(ptrdiff_t[]){0,0,p,p});
}

static const struct const_Matrix_R* constructor_grad_xyz_p
	(const int dir, const struct Multiarray_Operator*const grad_rst, const struct const_Multiarray_R*const metrics)
{
	assert(metrics->layout == 'C');
	const struct const_Matrix_R*const grad_r = grad_rst->data[0]->op_std;
	struct Matrix_R*const grad_xyz = constructor_zero_Matrix_R(grad_r->layout,grad_r->ext_0,grad_r->ext_1); // rtrnd.

	for (int d = 0 ; d < DIM; ++d) {
		const struct const_Vector_R metrics_V =
			{ .ext_0     = metrics->extents[0],
			  .owns_data = false,
			  .data = get_col_const_Multiarray_R(dir*DIM+d,metrics), };
		mm_diag_R('L',1.0,1.0,grad_rst->data[d]->op_std,&metrics_V,grad_xyz,false);
	}

	return (struct const_Matrix_R*) grad_xyz;
}

static void constructor_d_g_coef_f__d_s_coef_i
	(const int side_index, const int sol_index, const int ind_num_flux_2nd,
	 const struct const_Matrix_R*const jdet_n_fc, struct DG_Solver_Face_T*const dg_s_face, const bool collocated)
{
	const struct Face*const face            = (struct Face*) dg_s_face;
	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	assert(!face->boundary);

	const struct Operator* cv0_vs_fc_op = get_operator__cv0_vs_fc_T(sol_index,s_face);
	const struct const_Matrix_R* cv0_vs_fc = NULL;
	if (side_index == sol_index) {
		cv0_vs_fc = cv0_vs_fc_op->op_std;
	} else {
		const struct const_Vector_i* nc_fc = get_operator__nc_fc_T(sol_index,s_face);
		cv0_vs_fc = constructor_copy_permute_const_Matrix_R(cv0_vs_fc_op->op_std,nc_fc,'R'); // destructed
	}

	const struct const_Matrix_R*const m_inv_tw0_vt_fc =
		constructor_m_inv_tw0_vt_fc(side_index,dg_s_face,collocated); // destructed
	const struct const_Vector_R*const jn_fc_V = get_jn_fc_V(jdet_n_fc);

	struct Neigh_Info_DG*const neigh_info = &dg_s_face->neigh_info[side_index];

	Real sol_scale = get_sol_scale(side_index,sol_index,ind_num_flux_2nd);
	for (int d = 0; d < DIM; ++d) {
		const struct const_Matrix_R*const right =
			constructor_mm_diag_const_Matrix_R(sol_scale,cv0_vs_fc,&jn_fc_V[d],'L',false); // destructed

		destructor_const_Matrix_R(neigh_info->d_g_coef_f__d_s_coef[sol_index][d]);
		neigh_info->d_g_coef_f__d_s_coef[sol_index][d] =
			constructor_mm_const_Matrix_R('N','N',1.0,m_inv_tw0_vt_fc,right,'R'); // keep
		destructor_const_Matrix_R(right);
	}
	destructor_const_Matrix_R(m_inv_tw0_vt_fc);

	if (side_index != sol_index)
		destructor_const_Matrix_R(cv0_vs_fc);
}

static void constructor_d_g_coef_f__d_s_coef_b
	(const int ind_num_flux_2nd, const struct const_Matrix_R*const jdet_n_fc, struct DG_Solver_Face_T*const dg_s_face,
	 const struct Simulation*const sim)
{
	const int side_index = 0,
	          sol_index  = 0;

	const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) dg_s_face;

	const struct Operator*const tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);
	const struct Operator* cv0_vs_fc      = get_operator__cv0_vs_fc_T(sol_index,s_face);
	struct Neigh_Info_DG*const neigh_info = &dg_s_face->neigh_info[side_index];

	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	const ptrdiff_t ext_0 = tw0_vt_fc->op_std->ext_0,
	                ext_1 = cv0_vs_fc->op_std->ext_1;

	for (int d = 0; d < DIM; ++d) {
		destructor_const_Matrix_R(neigh_info->d_g_coef_f__d_s_coef[sol_index][d]);
		neigh_info->d_g_coef_f__d_s_coef[sol_index][d] =
			constructor_empty_const_Matrix_R('R',n_eq*ext_0,n_vr*ext_1); // keep
	}
	struct Matrix_R*const local_block = constructor_empty_Matrix_R('R',ext_0,ext_1); // destructed


	const bool collocated = sim->collocated;
	const struct const_Matrix_R*const m_inv_tw0_vt_fc =
		constructor_m_inv_tw0_vt_fc(side_index,dg_s_face,collocated); // destructed
	const struct const_Vector_R*const jn_fc_V = get_jn_fc_V(jdet_n_fc);

	const bool force_implicit = ( test_case->solver_method_curr == 'i' ? false : true );
	if (force_implicit)
		const_cast_c(&test_case->solver_method_curr,'i');

	struct Numerical_Flux_Input_T*const num_flux_i = constructor_Numerical_Flux_Input_T(sim); // destructed
	constructor_Numerical_Flux_Input_data_T(num_flux_i,s_face,sim); // destructed

	const struct Boundary_Value_T*const bv = &num_flux_i->bv_r;

	if (force_implicit)
		const_cast_c(&test_case->solver_method_curr,'e');

	for (int eq = 0; eq < n_eq; ++eq) {
	for (int vr = 0; vr < n_vr; ++vr) {

// See \ref numerical_flux_T.c for ind_ds_ds and usage.
UNUSED(bv);
UNUSED(jn_fc_V);
UNUSED(ind_num_flux_2nd);

// use get_sol_const to get scaling for u_i(nternal) and u_b(oundary) terms

		for (int d = 0; d < DIM; ++d) {
			const struct const_Matrix_R*const right = NULL;
EXIT_UNSUPPORTED;
//				constructor_mm_diag_const_Matrix_R(1.0,cv0_vs_fc->op_std,ds_ds_jn_fc_V,'L',false); // destructed
			mm_R('N','N',1.0,0.0,m_inv_tw0_vt_fc,right,local_block);
			destructor_const_Matrix_R(right);

			set_block_Matrix_R((struct Matrix_R*)neigh_info->d_g_coef_f__d_s_coef[sol_index][d],
			                   eq*ext_0,vr*ext_1,(struct const_Matrix_R*)local_block,0,0,ext_0,ext_1,'i');
		}
	}}
	destructor_Numerical_Flux_Input_data_T(num_flux_i);
	destructor_Numerical_Flux_Input_T(num_flux_i);

	destructor_Matrix_R(local_block);
	destructor_const_Matrix_R(m_inv_tw0_vt_fc);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Matrix_R* constructor_m_inv_tw0_vt_fc
	(const int side_index, struct DG_Solver_Face_T*const dg_s_face, const bool collocated)
{
	const struct Face*const face                   = (struct Face*) dg_s_face;
	const struct Solver_Face_T*const s_face        = (struct Solver_Face_T*) dg_s_face;
	const struct Solver_Volume_T*const s_vol       = (struct Solver_Volume_T*) face->neigh_info[side_index].volume;
	const struct DG_Solver_Volume_T*const dg_s_vol = (struct DG_Solver_Volume_T*) s_vol;

	const struct Operator*const tw0_vt_fc = get_operator__tw0_vt_fc_T(side_index,s_face);

	const struct const_Matrix_R* m_inv_tw = NULL;
	if (!collocated) {
		m_inv_tw = constructor_mm_const_Matrix_R('N','N',1.0,dg_s_vol->m_inv,tw0_vt_fc->op_std,'R'); // returned
	} else {
		const struct const_Vector_R*const w_vc  = get_operator__w_vc__s_e_T(s_vol);
		const struct const_Vector_R jdet_vc     = interpret_const_Multiarray_as_Vector_R(s_vol->jacobian_det_vc);
		const struct const_Vector_R*const wJ_vc = constructor_dot_mult_const_Vector_R(1.0,w_vc,&jdet_vc,1); // dest.
		m_inv_tw = constructor_mm_diag_const_Matrix_R(1.0,tw0_vt_fc->op_std,wJ_vc,'L',true); // returned
		destructor_const_Vector_R(wJ_vc);
	}
	return m_inv_tw;
}

static struct const_Vector_R* get_jn_fc_V (const struct const_Matrix_R*const jdet_n_fc)
{
	static struct Vector_R jn_fc_V[DIM];

	for (int d = 0; d < DIM; ++d) {
		jn_fc_V[d].ext_0 = jdet_n_fc->ext_0;
		jn_fc_V[d].data  = (double*)get_col_const_Matrix_R(d,jdet_n_fc);
	}

	return (struct const_Vector_R*) jn_fc_V;
}

static Real get_sol_scale (const int side_index, const int sol_index, const int ind_num_flux_2nd)
{
	switch (ind_num_flux_2nd) {
	case NUM_FLUX_BR2_STABLE: // fallthrough
	case NUM_FLUX_CDG2:
		// Note: Central numerical solution (u_num == 0.5*(u_L + u_R)).
		return 0.5*( side_index == sol_index ? -1.0 : 1.0 );
		break;
	default:
		EXIT_ERROR("Unsupported: %d",ind_num_flux_2nd);
		break;
	}
}
