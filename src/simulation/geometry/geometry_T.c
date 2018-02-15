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
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_mesh.h"
#include "definitions_intrusive.h"
#include "definitions_visualization.h"


#include "def_templates_geometry.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"

#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

#define OUTPUT_GEOMETRY false ///< Flag for whether the geometry should be output for visualization.

/** \brief Pointer to functions computing \ref Solver_Volume_T::geom_coef.
 *  \param sim   \ref Simulation.
 *  \param s_vol \ref Solver_Volume_T.
 */
typedef void (*compute_geom_coef_fptr_T)
	(const struct Simulation*const sim,
	 struct Solver_Volume_T*const s_vol
	);

/// \brief Compute the face normal vectors at the nodes corresponding to the given face metrics.
static void compute_normals_T
	(const int ind_lf,                             ///< Defined for \ref compute_unit_normals_and_det_T.
	 const struct const_Multiarray_R* normals_ref, ///< Defined for \ref compute_unit_normals_and_det_T.
	 const struct const_Multiarray_R* metrics_f,   ///< Defined for \ref compute_unit_normals_and_det_T.
	 struct Multiarray_R* normals_f                ///< Defined for \ref compute_unit_normals_and_det_T.
	);

/** \brief Set the appropriate function pointer for computing \ref Solver_Volume_T::geom_coef.
 *  \return See brief. */
static compute_geom_coef_fptr_T set_fptr_geom_coef_T
	(const int domain_type,   ///< \ref Simulation::domain_type.
	 const bool volume_curved ///< \ref Volume::curved.
	);

/** \brief See return.
 *  \return The permutation required for conversion to the standard Jacobian ordering from the transposed ordering. */
static const ptrdiff_t* set_jacobian_permutation
	(const int d ///< The dimension
	);

/// \brief Compute the determinant of the geometry mapping Jacobian.
static void compute_detJV_T
	(struct const_Multiarray_R* jacobian, ///< Multiarray containing the Jacobian terms
	 struct Multiarray_R* jacobian_det    ///< Multiarray to be set to contain the determinant.
	);

/// \brief Compute the cofactor matrix entries of the geometry mapping Jacobian (referred to as the metrics).
static void compute_cofactors_T
	(struct const_Multiarray_R* jacobian, ///< Multiarray containing the Jacobian terms
	 struct Multiarray_R* metrics         ///< Multiarray to be set to contain the metric terms.
	);

/** \brief Constructor for the "xyz" coordinates at the 'f'ace 'c'ubature nodes with a possible correction ensuring that
 *         they are placed on the exact domain boundary in the case of curved boundary faces.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_fc_with_exact_boundary
	(const struct Solver_Face_T*const s_face, ///< \ref Solver_Face_T.
	 const struct Simulation*const sim        ///< \ref Simulation.
	);

/** \brief Compute the face unit normal vectors at the nodes corresponding to the given face metrics.
 *  The l^2 norm of the initially un-normalized normal vector at each of the nodes is stored in `jacobian_det`. This is
 *  done in accordance with the definition (see (eq. (B.6), \cite Zwanenburg2016)).
 */
static void compute_unit_normals_and_det_T
	(const int ind_lf,                             ///< \ref Face::Neigh_Info::ind_lf in \ref Face.
	 const struct const_Multiarray_R* normals_ref, ///< \ref Element::normals.
	 const struct const_Multiarray_R* metrics_f,   /**< \ref Solver_Volume_T::metrics_vm interpolated to the face
	                                                *   nodes. */
	 struct Multiarray_R* normals_f,               ///< \ref Multiarray_T\* in which to store the face normals.
	 struct Multiarray_R* jacobian_det_f           /**< \ref Multiarray_T\* in which to store the face jacobian
	                                                *   determinants. */
	);

// Interface functions ********************************************************************************************** //

void set_up_solver_geometry_T (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		compute_geometry_volume_T((struct Solver_Volume_T*)curr,sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		compute_geometry_face_T((struct Solver_Face_T*)curr,sim);

#if TYPE_RC == TYPE_REAL
	if (OUTPUT_GEOMETRY) {
		output_visualization(sim,VIS_GEOM_VOLUMES);
		EXIT_UNSUPPORTED;
	}
#endif
}

void compute_unit_normals_T
	(const int ind_lf, const struct const_Multiarray_R* normals_ref, const struct const_Multiarray_R* metrics_f,
	 struct Multiarray_R* normals_f)
{
	compute_normals_T(ind_lf,normals_ref,metrics_f,normals_f);
	normalize_Multiarray_R(normals_f,"L2",false,NULL);
}

void compute_geometry_volume_T (struct Solver_Volume_T* s_vol, const struct Simulation *sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	compute_geom_coef_fptr_T compute_geom_coef = set_fptr_geom_coef_T(sim->domain_type,vol->curved);
	compute_geom_coef(sim,s_vol);

	const int d = ((struct const_Element*)g_e)->d;

	const struct const_Multiarray_R*const geom_coef = s_vol->geom_coef;

	const int p = s_vol->p_ref;

	const bool curved = vol->curved;
	const int p_g = ( curved ? p : 1 );

	struct Ops {
		const struct Multiarray_Operator cv1_vg_vm;
		const struct Multiarray_Operator cv1_vg_vc;
		const struct Operator* vv0_vm_vc;
	} ops = { .cv1_vg_vm = set_MO_from_MO(g_e->cv1_vg_vm[curved],1,(ptrdiff_t[]){0,0,p_g,p_g}),
	          .cv1_vg_vc = set_MO_from_MO(g_e->cv1_vg_vc[curved],1,(ptrdiff_t[]){0,0,p,p_g}),
	          .vv0_vm_vc = get_Multiarray_Operator(g_e->vv0_vm_vc[curved],(ptrdiff_t[]){0,0,p,p_g}), };

	const ptrdiff_t n_vm = ops.cv1_vg_vm.data[0]->op_std->ext_0,
	                n_vc = ops.cv1_vg_vc.data[0]->op_std->ext_0;

	struct Multiarray_R* jacobian_vm = constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){n_vm,d,d}), // destructed
	                   * jacobian_vc = constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){n_vc,d,d}); // destructed

	for (ptrdiff_t row = 0; row < d; ++row) {
		mm_NN1C_Operator_Multiarray_R(ops.cv1_vg_vm.data[row],geom_coef,jacobian_vm,op_format,2,NULL,&row);
		mm_NN1C_Operator_Multiarray_R(ops.cv1_vg_vc.data[row],geom_coef,jacobian_vc,op_format,2,NULL,&row);
	}

	const ptrdiff_t* perm = set_jacobian_permutation(d);
	permute_Multiarray_R(jacobian_vc,perm,jacobian_vc->layout);

	compute_detJV_T((struct const_Multiarray_R*)jacobian_vc,(struct Multiarray_R*)s_vol->jacobian_det_vc);
	compute_cofactors_T((struct const_Multiarray_R*)jacobian_vm,(struct Multiarray_R*)s_vol->metrics_vm);

	destructor_Multiarray_R(jacobian_vm);
	destructor_Multiarray_R(jacobian_vc);

	const struct const_Multiarray_R* met_vm = s_vol->metrics_vm;

	resize_Multiarray_R((struct Multiarray_R*)s_vol->metrics_vc,3,(ptrdiff_t[]){n_vc,d,d});
	mm_NN1C_Operator_Multiarray_R(
		ops.vv0_vm_vc,met_vm,(struct Multiarray_R*)s_vol->metrics_vc,op_format,met_vm->order,NULL,NULL);
}

void compute_geometry_face_T (struct Solver_Face_T* s_face, struct Simulation* sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	assert((s_face->cub_type == 's') || (s_face->cub_type == 'c'));

	struct Face* face             = (struct Face*) s_face;
	struct Volume* vol            = face->neigh_info[0].volume;
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	struct const_Element* e            = (struct const_Element*) g_e;

	struct Ops {
		const struct Operator* cv0_vg_fc;
		const struct Operator* vv0_vm_fc;
	} ops = { .cv0_vg_fc = NULL,
	          .vv0_vm_fc = NULL, };

	const int ind_lf = face->neigh_info[0].ind_lf;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved_f = (s_face->cub_type == 's' ? 0 : 1);
	if (!vol->curved) {
		ops.cv0_vg_fc = get_Multiarray_Operator(g_e->cv0_vgs_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,1});
		ops.vv0_vm_fc = get_Multiarray_Operator(g_e->vv0_vms_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,1});
	} else {
		ops.cv0_vg_fc = get_Multiarray_Operator(g_e->cv0_vgc_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,p_v});
		ops.vv0_vm_fc = get_Multiarray_Operator(g_e->vv0_vmc_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,p_v});
	}

	destructor_const_Multiarray_R(s_face->xyz_fc);
	const_constructor_move_const_Multiarray_R(&s_face->xyz_fc,
	                                          constructor_xyz_fc_with_exact_boundary(s_face,sim)); // keep

	const struct const_Multiarray_R* m_vm = s_vol->metrics_vm;
	const struct const_Multiarray_R* metrics_fc =
		constructor_mm_NN1_Operator_const_Multiarray_R(ops.vv0_vm_fc,m_vm,'C',op_format,m_vm->order,NULL); // destructed

	compute_unit_normals_and_det_T(ind_lf,e->normals,metrics_fc,
		(struct Multiarray_R*)s_face->normals_fc,(struct Multiarray_R*)s_face->jacobian_det_fc);

	destructor_const_Multiarray_R(metrics_fc);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref compute_geom_coef_fptr_T for straight volumes.
static void compute_geom_coef_straight_T
	(const struct Simulation*const sim, ///< See brief.
	 struct Solver_Volume_T*const s_vol ///< See brief.
	);

/// \brief Version of \ref compute_geom_coef_fptr_T for curved volumes using blending.
static void compute_geom_coef_blended_T
	(const struct Simulation*const sim, ///< See brief.
	 struct Solver_Volume_T*const s_vol ///< See brief.
	);

/// \brief Version of \ref compute_geom_coef_fptr_T for curved volumes using the parametric mapping.
static void compute_geom_coef_parametric_T
	(const struct Simulation*const sim, ///< See brief.
	 struct Solver_Volume_T*const s_vol ///< See brief.
	);

static void compute_normals_T
	(const int ind_lf, const struct const_Multiarray_R* normals_ref, const struct const_Multiarray_R* metrics_f,
	 struct Multiarray_R* normals_f)
{
	const int order_n = 2,
	          order_m = 3;

	const ptrdiff_t* exts_m = metrics_f->extents;

	assert(normals_ref->order == order_n);
	assert(metrics_f->order == order_m);
	assert(exts_m[1] == exts_m[2]);
	assert(normals_f->order == order_n);
	assert(metrics_f->layout == 'C');

	const ptrdiff_t n_vals = exts_m[0],
	                d      = exts_m[1];

	resize_Multiarray_R(normals_f,order_n,(ptrdiff_t[]){n_vals,d});
	set_to_value_Multiarray_R(normals_f,0.0);

	Real* normals_d = normals_f->data;
	const Real* normal_ref = get_row_const_Multiarray_R(ind_lf,normals_ref);
	const Real* metrics_d  = metrics_f->data;
	for (ptrdiff_t dim_0 = 0; dim_0 < d; ++dim_0) {
	for (ptrdiff_t n = 0; n < n_vals; ++n) {
		for (ptrdiff_t dim_1 = 0; dim_1 < d; ++dim_1)
			*normals_d += normal_ref[dim_1]*metrics_d[n_vals*(dim_0+d*dim_1)+n];
		++normals_d;
	}}
	normals_f->layout = 'C';

	transpose_Multiarray_R(normals_f,true);
}

static compute_geom_coef_fptr_T set_fptr_geom_coef_T (const int domain_type, const bool volume_curved)
{
	if (domain_type == DOM_STRAIGHT) {
		return compute_geom_coef_straight_T;
	} else if (domain_type == DOM_BLENDED) {
		if (!volume_curved)
			return compute_geom_coef_straight_T;
		else
			return compute_geom_coef_blended_T;
	} else if (domain_type == DOM_PARAMETRIC) {
		return compute_geom_coef_parametric_T;
	}

	EXIT_ERROR("Unsupported domain_type: %d\n",domain_type);
}

static const ptrdiff_t* set_jacobian_permutation (const int d)
{
	switch (d) {
	case 1:
		return NULL;
		break;
	case 2: {
		static const ptrdiff_t perm_2d[] = { 0, 2, 1, 3, };
		return perm_2d;
		break;
	} case 3: {
		static const ptrdiff_t perm_3d[] = { 0, 3, 6, 1, 4, 7, 2, 5, 8, };
		return perm_3d;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",d);
		break;
	}
	return NULL;
}

static void compute_detJV_T (struct const_Multiarray_R* jacobian, struct Multiarray_R* jacobian_det)
{
	const int order_j  = 3,
	          order_dj = 1;

	const ptrdiff_t* exts_j = jacobian->extents;

	assert(jacobian_det->order == order_dj);
	assert(jacobian->order == order_j);
	assert(exts_j[1] == exts_j[2]);

	const ptrdiff_t n_vals = exts_j[0],
	                d      = exts_j[1];

	resize_Multiarray_R(jacobian_det,order_dj,(ptrdiff_t[]){n_vals});
	Real* j_det = jacobian_det->data;

	switch (d) {
	case 1: {
		const Real* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})];
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			j_det[i] = x_r[i];
		break;
	} case 2: {
		const Real* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
		          * x_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,1})],
		          * y_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,0})],
		          * y_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,1})];
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			j_det[i] = x_r[i]*y_s[i]-x_s[i]*y_r[i];
		break;
	} case 3: {
		const Real* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
		          * x_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,1})],
		          * x_t = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,2})],
		          * y_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,0})],
		          * y_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,1})],
		          * y_t = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,2})],
		          * z_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){2,0})],
		          * z_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){2,1})],
		          * z_t = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){2,2})];
		for (ptrdiff_t i = 0; i < n_vals; ++i) {
			j_det[i] = x_r[i]*(y_s[i]*z_t[i]-y_t[i]*z_s[i])
			          -x_s[i]*(y_r[i]*z_t[i]-y_t[i]*z_r[i])
			          +x_t[i]*(y_r[i]*z_s[i]-y_s[i]*z_r[i]);
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %td\n",d);
		break;
	}

	for (ptrdiff_t i = 0; i < n_vals; ++i)
		assert(j_det[i] > 0.0);
}

static void compute_cofactors_T (struct const_Multiarray_R* jacobian, struct Multiarray_R* metrics)
{
	const int order_j = 3,
	          order_m = 3;

	const ptrdiff_t* exts_j = jacobian->extents,
	               * exts_m = metrics->extents;

	assert(metrics->order == order_m);
	assert(jacobian->order == order_j);
	assert(exts_j[1] == exts_j[2]);

	const ptrdiff_t n_vals = exts_j[0],
	                d      = exts_j[1];

	resize_Multiarray_R(metrics,order_m,(ptrdiff_t[]){n_vals,d,d});
	switch (d) {
	case 1:
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			metrics->data[i] = 1.0;
		break;
	case 2: {
		const Real* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
		          * x_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,1})],
		          * y_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,0})],
		          * y_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,1})];
		Real* m_00 = &metrics->data[compute_index_sub_container(order_m,1,exts_m,(ptrdiff_t[]){0,0})],
		    * m_01 = &metrics->data[compute_index_sub_container(order_m,1,exts_m,(ptrdiff_t[]){0,1})],
		    * m_10 = &metrics->data[compute_index_sub_container(order_m,1,exts_m,(ptrdiff_t[]){1,0})],
		    * m_11 = &metrics->data[compute_index_sub_container(order_m,1,exts_m,(ptrdiff_t[]){1,1})];
		for (ptrdiff_t i = 0; i < n_vals; ++i) {
			m_00[i] =  y_s[i];
			m_01[i] = -y_r[i];
			m_10[i] = -x_s[i];
			m_11[i] =  x_r[i];
		}
		break;
	} case 3: {
		EXIT_ADD_SUPPORT; // curl-form. See comments in \ref compute_geometry_volume_T_T.
		break;
	} default:
		EXIT_ERROR("Unsupported: %td\n",d);
		break;
	}
}

static const struct const_Multiarray_R* constructor_xyz_fc_with_exact_boundary
	(const struct Solver_Face_T*const s_face, const struct Simulation*const sim)
{
	struct Face*const face             = (struct Face*) s_face;
	struct Volume*const vol            = face->neigh_info[0].volume;
	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	const struct Operator* cv0_vg_fc = NULL;

	const int ind_lf = face->neigh_info[0].ind_lf;
	const int p_v = s_vol->p_ref,
	          p_f = s_face->p_ref;

	const int curved_f = (s_face->cub_type == 's' ? 0 : 1);
	if (!vol->curved)
		cv0_vg_fc = get_Multiarray_Operator(g_e->cv0_vgs_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,1});
	else
		cv0_vg_fc = get_Multiarray_Operator(g_e->cv0_vgc_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,p_v});

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_f = 'd';

	const struct const_Multiarray_R*const g_coef = s_vol->geom_coef;
	struct Multiarray_R*const xyz_fc = constructor_mm_NN1_Operator_Multiarray_R
		(cv0_vg_fc,(struct Multiarray_R*)g_coef,'C',op_f,g_coef->order,NULL); // returned

	if (is_face_bc_curved(face->bc)) {
		assert(face->neigh_info[0].ind_href == 0);

		const char ce_type = 'f',
		           n_type  = 'c';
		const int p_g = s_vol->p_ref;

		struct Boundary_Comp_Elem_Data_T b_ce_d =
			constructor_static_Boundary_Comp_Elem_Data_T(ce_type,n_type,p_g,s_vol); // destructed
		set_Boundary_Comp_Elem_operators_T(&b_ce_d,s_vol,ce_type,n_type,p_g,ind_lf);

		struct Matrix_R xyz_fc_M = interpret_Multiarray_as_Matrix_R(xyz_fc);
		const struct const_Matrix_R*const xyz_fc_diff =
			constructor_xyz_surf_diff_T(&b_ce_d,(struct const_Matrix_R*)&xyz_fc_M,s_vol,n_type,sim); // d
		add_in_place_Matrix_R(1.0,&xyz_fc_M,xyz_fc_diff);

		destructor_const_Matrix_R(xyz_fc_diff);
		destructor_static_Boundary_Comp_Elem_Data_T(&b_ce_d);
	}

	return (const struct const_Multiarray_R*) xyz_fc;
}

static void compute_unit_normals_and_det_T
	(const int ind_lf, const struct const_Multiarray_R* normals_ref, const struct const_Multiarray_R* metrics_f,
	 struct Multiarray_R* normals_f, struct Multiarray_R* jacobian_det_f)
{
	compute_normals_T(ind_lf,normals_ref,metrics_f,normals_f);
	normalize_Multiarray_R(normals_f,"L2",true,jacobian_det_f);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the high-order straight geometry values.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_s_ho
	(const struct Solver_Volume_T*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

/** \brief Constructor for the high-order geometry coefficients from the input values.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_geom_coef_ho
	(const struct const_Multiarray_R* geom_val, ///< The geometry values.
	 const struct Solver_Volume_T*const s_vol,  ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim          ///< \ref Simulation.
	);

static void compute_geom_coef_straight_T (const struct Simulation*const sim, struct Solver_Volume_T*const s_vol)
{
	struct Volume* vol = (struct Volume*) s_vol;
	destructor_const_Multiarray_R(s_vol->geom_coef);

	if (strstr(sim->basis_geom,"lagrange") || strstr(sim->basis_geom,"bezier")) {
		const_constructor_copy_Multiarray_R(&s_vol->geom_coef,vol->xyz_ve); // keep
		if (s_vol->geom_coef->layout != 'C')
			transpose_Multiarray_R((struct Multiarray_R*)s_vol->geom_coef,true);
	} else if (strstr(sim->basis_geom,"nurbs")) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_ERROR("Unsupported: '%s'.",sim->basis_geom);
	}
}

static void compute_geom_coef_blended_T (const struct Simulation*const sim, struct Solver_Volume_T*const s_vol)
{
	const struct const_Multiarray_R* xyz_s = constructor_xyz_s_ho(s_vol,sim); // destructed

	const struct const_Multiarray_R* xyz = constructor_xyz_blended_T('g',xyz_s,s_vol,sim); // destructed
	destructor_const_Multiarray_R(xyz_s);

	destructor_const_Multiarray_R(s_vol->geom_coef);
	const_constructor_move_const_Multiarray_R(&s_vol->geom_coef,constructor_geom_coef_ho(xyz,s_vol,sim)); // keep
	destructor_const_Multiarray_R(xyz);
}

static void compute_geom_coef_parametric_T (const struct Simulation*const sim, struct Solver_Volume_T*const s_vol)
{
	const struct const_Multiarray_R* xyz_s = constructor_xyz_s_ho(s_vol,sim); // destructed

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const struct const_Multiarray_R* xyz = test_case->constructor_xyz(0,xyz_s,s_vol,sim); // destructed
	destructor_const_Multiarray_R(xyz_s);

/// \todo add setup_bezier_mesh: Likely make this a separate function.

	destructor_const_Multiarray_R(s_vol->geom_coef);
	const_constructor_move_const_Multiarray_R(&s_vol->geom_coef,constructor_geom_coef_ho(xyz,s_vol,sim)); // keep
	destructor_const_Multiarray_R(xyz);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Multiarray_R* constructor_xyz_s_ho
	(const struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	const bool curved = vol->curved;
	const int p       = s_vol->p_ref;
	assert(curved == true);

	const struct Operator* vv0_vv_vg = get_Multiarray_Operator(g_e->vv0_vv_vg[curved],(ptrdiff_t[]){0,0,p,1});
	const struct const_Multiarray_R* xyz_ve = vol->xyz_ve;

	return constructor_mm_NN1_Operator_const_Multiarray_R(vv0_vv_vg,xyz_ve,'C',op_format,xyz_ve->order,NULL);
}

static const struct const_Multiarray_R* constructor_geom_coef_ho
	(const struct const_Multiarray_R* geom_val, const struct Solver_Volume_T*const s_vol,
	 const struct Simulation*const sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int p = s_vol->p_ref;

	const struct Operator* vc0_vgc_vgc = get_Multiarray_Operator(g_e->vc0_vgc_vgc,(ptrdiff_t[]){0,0,p,p});

	return constructor_mm_NN1_Operator_const_Multiarray_R(vc0_vgc_vgc,geom_val,'C',op_format,geom_val->order,NULL);
}
