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
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dg.h"

#include "def_templates_compute_face_rlhs.h"
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

/** \brief Compute the scaling associated with the \ref Solver_Volume_T::sol_coef contribution to the current face weak
 *         gradient term in the integrated by parts twice form. */
static Real get_sol_scale
	(const int side_index,      ///< The side index of the face.
	 const int sol_index,       ///< The side index of the solution.
	 const int ind_num_flux_2nd ///< The second component of \ref Test_Case_T::ind_num_flux.
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
				const struct const_Matrix_R*const d_ibp2 =
					constructor_mm_const_Matrix_R('N','N',1.0,tw0_vt_vc->op_std,grad_xyz,'R'); // destructed
				d_g_coef_v__d_s_coef =
					constructor_mm_const_Matrix_R('N','N',1.0,dg_s_vol->m_inv,d_ibp2,'R'); // keep
				destructor_const_Matrix_R(d_ibp2);
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
//		struct DG_Solver_Face_T*const dg_s_face = (struct DG_Solver_Face_T*) curr;

		const struct Operator* tw0_vt_fc[2] = { get_operator__tw0_vt_fc_T(0,s_face), NULL, };
		const struct Operator* cv0_vs_fc[2] = { get_operator__cv0_vs_fc_T(0,s_face), NULL, };

		const struct const_Vector_R jdet_fc    = interpret_const_Multiarray_as_Vector_R(s_face->jacobian_det_fc);
		const struct const_Matrix_R normals_fc = interpret_const_Multiarray_as_Matrix_R(s_face->normals_fc);

		struct Matrix_R*const jdet_normals_fc =
			constructor_mm_diag_Matrix_R(1.0,&normals_fc,&jdet_fc,'L',false); // destructed

		const struct Test_Case_T*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
		Real sol_scale = get_sol_scale(0,0,test_case->ind_num_flux[1]);
UNUSED(sol_scale);

		if (!face->boundary) {
			tw0_vt_fc[1] = get_operator__tw0_vt_fc_T(1,s_face);
			cv0_vs_fc[1] = get_operator__cv0_vs_fc_T(1,s_face);

			// Note: The indices here are the indices of the destination face.
			const struct const_Vector_i* nc_fc[2] = { get_operator__nc_fc_T(1,s_face),
			                                          get_operator__nc_fc_T(0,s_face), };

			const struct const_Matrix_R*const cv0_vs_fc_M =
				constructor_copy_permute_const_Matrix_R(cv0_vs_fc[0]->op_std,nc_fc[0],'R'); // destructed
UNUSED(cv0_vs_fc_M);
		}

		// function to compute scaling value associated with the numerical-interpolated flux (dependent upon
		// NUM_FLUX.
		destructor_Matrix_R(jdet_normals_fc);
EXIT_UNSUPPORTED;
	}
}

// Level 1 ********************************************************************************************************** //

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

static Real get_sol_scale (const int side_index, const int sol_index, const int ind_num_flux_2nd)
{
	switch (ind_num_flux_2nd) {
	case NUM_FLUX_BR2_STABLE: // fallthrough
	case NUM_FLUX_CDG2:
		// Note: Central numerical solution.
		return 0.5*( side_index == sol_index ? -1.0 : 1.0 );
		break;
	default:
		EXIT_ERROR("Unsupported: %d",ind_num_flux_2nd);
		break;
	}
}
