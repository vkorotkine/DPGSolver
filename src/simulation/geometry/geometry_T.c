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

#include "def_templates_geometry_normals.h"
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

/// \brief Version of \ref compute_geom_coef_fptr_T for degree 1 coefficients.
static void compute_geom_coef_p1
	(const struct Simulation*const sim, ///< See brief.
	 struct Solver_Volume_T*const s_vol ///< See brief.
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

/** \brief Constructor for the "xyz" coordinates at the 'f'ace 'c'ubature nodes interpolated from the volume geometry
 *         nodes.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_fc
	(const struct Solver_Face_T*const s_face, ///< \ref Solver_Face_T.
	 const struct Simulation*const sim        ///< \ref Simulation.
	);

/** \brief Constructor for the "xyz" coordinates at the 'f'ace 'c'ubature nodes with a possible correction ensuring that
 *         they are placed on the exact domain boundary in the case of curved boundary faces.
 *  \return See brief.
 *
 *  If the correction is not used, xyz-coordinate dependent boundary conditions on curved boundary faces are set based
 *  on the approximate geometry. For Poisson solutions, this results in optimal convergence even when isoparametric
 *  geometry is not used. For Advection solutions however, this correction seems to be introducing an error causing the
 *  loss of free-stream preservation when a non-constant advection velocity is used.
 *
 *  \todo Investigate further and update comments above.
 */
static const struct const_Multiarray_R* constructor_xyz_fc_on_exact_boundary
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

/// \brief Compute \ref Solver_Face_T::vol_jacobian_det_fc.
static void compute_vol_jacobian_det_fc_T
	(struct Solver_Face_T*const s_face ///< Standard.
	);

/** \brief Set data relating to the p1 geometry of the input \ref Solver_Volume_T.
 *
 *  This function sets:
 *  - \ref Solver_Volume_T::metrics_vm_p1;
 */
static void compute_geometry_volume_p1_T
	(struct Solver_Volume_T*const s_vol ///< Standard.
	);

/** \brief Set data relating to the p1 geometry of the input \ref Solver_Face_T.
 *
 *  This function sets:
 *  - \ref Solver_Face_T::normals_p1;
 *  - \ref Solver_Face_T::jacobian_det_p1;
 */
static void compute_geometry_face_p1_T
	(struct Solver_Face_T*const s_face ///< Standard.
	);

// Interface functions ********************************************************************************************** //

void set_up_solver_geometry_T (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		compute_geometry_volume_T(true,(struct Solver_Volume_T*)curr,sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		compute_geometry_face_T((struct Solver_Face_T*)curr,sim);

	correct_for_exact_normals_T(sim);

#if TYPE_RC == TYPE_REAL
	if (OUTPUT_GEOMETRY) {
		output_visualization(sim,VIS_GEOM_VOLUMES);
		EXIT_ERROR("Disable output to continue.\n");
	}
#endif
}

void set_up_solver_geometry_p1_T (struct Simulation*const sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		compute_geometry_volume_p1_T((struct Solver_Volume_T*)curr);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		compute_geometry_face_p1_T((struct Solver_Face_T*)curr);
}

void compute_unit_normals_T
	(const int ind_lf, const struct const_Multiarray_R* normals_ref, const struct const_Multiarray_R* metrics_f,
	 struct Multiarray_R* normals_f)
{
	compute_normals_T(ind_lf,normals_ref,metrics_f,normals_f);
	normalize_Multiarray_R(normals_f,"L2",false,NULL);
}

void compute_geometry_volume_T
	(const bool recompute_geom_coef, struct Solver_Volume_T* s_vol, const struct Simulation *sim)
{
	const char op_format = get_set_op_format(0);

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	if (recompute_geom_coef) {
		compute_geom_coef_fptr_T compute_geom_coef = set_fptr_geom_coef_T(sim->domain_type,vol->curved);
		compute_geom_coef_p1(sim,s_vol);
		compute_geom_coef(sim,s_vol);
	}

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

void compute_geometry_face_T (struct Solver_Face_T* s_face, const struct Simulation*const sim)
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
	destructor_const_Multiarray_R(s_face->xyz_fc_ex_b);
	const_constructor_move_const_Multiarray_R(&s_face->xyz_fc,constructor_xyz_fc(s_face,sim)); // keep
	const_constructor_move_const_Multiarray_R(&s_face->xyz_fc_ex_b,
	                                          constructor_xyz_fc_on_exact_boundary(s_face,sim)); // keep

	const struct const_Multiarray_R* m_vm = s_vol->metrics_vm;
	const struct const_Multiarray_R* metrics_fc =
		constructor_mm_NN1_Operator_const_Multiarray_R(ops.vv0_vm_fc,m_vm,'C',op_format,m_vm->order,NULL); // d.

	compute_unit_normals_and_det_T(ind_lf,e->normals,metrics_fc,
		(struct Multiarray_R*)s_face->normals_fc,(struct Multiarray_R*)s_face->jacobian_det_fc);

	compute_vol_jacobian_det_fc_T(s_face);

	destructor_const_Multiarray_R(metrics_fc);
}

const struct const_Multiarray_R* constructor_xyz_s_ho_T
	(const char ve_rep, const struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
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
	assert(ve_rep == 'v' || ve_rep == 'g');
	const struct const_Multiarray_R* xyz_ve = NULL;
	if (ve_rep == 'v') {
		xyz_ve = vol->xyz_ve;
	} else if (ve_rep == 'g') {
		const int pg = ( curved ? p : 1 );
		const struct Operator* cv0_vg_vv = get_Multiarray_Operator(g_e->cv0_vg_vv[curved],(ptrdiff_t[]){0,0,1,pg});
		const struct const_Multiarray_R*const g_coef = s_vol->geom_coef;
		xyz_ve = constructor_mm_NN1_Operator_const_Multiarray_R(
			cv0_vg_vv,g_coef,'C',op_format,g_coef->order,NULL); // destructed
	}

	const struct const_Multiarray_R*const xyz_s =
		constructor_mm_NN1_Operator_const_Multiarray_R(vv0_vv_vg,xyz_ve,'C',op_format,xyz_ve->order,NULL);

	if (ve_rep == 'g')
		destructor_const_Multiarray_R(xyz_ve);

	return xyz_s;
}

const struct const_Multiarray_R* constructor_geom_coef_ho_T
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

void correct_for_exact_normals_T (const struct Simulation*const sim)
{
	if (!using_exact_normals() && !using_exact_normals_for_boundary())
		return;

	correct_for_exact_normal_fptr_T correct_for_exact_normal = set_correct_for_exact_normal_fptr_T(sim);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		correct_for_exact_normal((struct Solver_Face_T*)curr);
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

	for (ptrdiff_t i = 0; i < n_vals; ++i) {
#ifndef NDEBUG
		if (j_det[i] < 0.0)
			print_Multiarray_R(jacobian_det);
#endif
		assert(j_det[i] > 0.0);
	}
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

static const struct const_Multiarray_R* constructor_xyz_fc
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
	return constructor_mm_NN1_Operator_const_Multiarray_R(cv0_vg_fc,g_coef,'C',op_f,g_coef->order,NULL);
}

static const struct const_Multiarray_R* constructor_xyz_fc_on_exact_boundary
	(const struct Solver_Face_T*const s_face, const struct Simulation*const sim)
{
	struct Face*const face             = (struct Face*) s_face;
	struct Volume*const vol            = face->neigh_info[0].volume;
	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;

	const int ind_lf = face->neigh_info[0].ind_lf;

	struct Multiarray_R*const xyz_fc = constructor_copy_Multiarray_R((struct Multiarray_R*)s_face->xyz_fc);
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
			constructor_xyz_surf_diff_T(&b_ce_d,(struct const_Matrix_R*)&xyz_fc_M,s_vol,n_type,false,sim); // d
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

static void compute_vol_jacobian_det_fc_T (struct Solver_Face_T*const s_face)
{
	const char op_format = 'd'; // Tensor-product operators potentially not provided.

	struct Face*const face                   = (struct Face*) s_face;
	const struct Volume*const vol            = (struct Volume*) face->neigh_info[0].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) face->neigh_info[0].volume;
	const struct Geometry_Element*const g_e  = &((struct Solver_Element*)vol->element)->g_e;

	const bool curved = vol->curved;
	const int p_g = ( curved ? s_vol->p_ref : 1 ),
	          p = s_face->p_ref,
	          ind_lf = face->neigh_info[0].ind_lf;

	const struct Multiarray_Operator* g_e__cv1_vg_fc = (!curved ? g_e->cv1_vgs_fc[curved] : g_e->cv1_vgc_fc[curved]);
//print_Multiarray_Operator(g_e__cv1_vg_fc);
printf("%d %p\n",curved,g_e__cv1_vg_fc);
printf("%d %d %d% d %d\n",ind_lf,0,0,p,p_g);
	const struct Multiarray_Operator cv1_vg_fc = set_MO_from_MO(g_e__cv1_vg_fc,1,(ptrdiff_t[]){ind_lf,0,0,p,p_g});

	const ptrdiff_t n_fc = cv1_vg_fc.data[0]->op_std->ext_0;

	struct Multiarray_R*const jac_fc = constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){n_fc,DIM,DIM}); // dest.

	for (ptrdiff_t row = 0; row < DIM; ++row)
		mm_NN1C_Operator_Multiarray_R(cv1_vg_fc.data[row],s_vol->geom_coef,jac_fc,op_format,2,NULL,&row);

	const ptrdiff_t*const perm = set_jacobian_permutation(DIM);
	permute_Multiarray_R(jac_fc,perm,jac_fc->layout);
	compute_detJV_T((struct const_Multiarray_R*)jac_fc,(struct Multiarray_R*)s_face->vol_jacobian_det_fc);
}

static void compute_geometry_volume_p1_T (struct Solver_Volume_T*const s_vol)
{
	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	const int d = ((struct const_Element*)g_e)->d;

	const struct Multiarray_Operator cv1_vg_vm = set_MO_from_MO(g_e->cv1_vg_vm[0],1,(ptrdiff_t[]){0,0,1,1});

	const ptrdiff_t n_vm = cv1_vg_vm.data[0]->op_std->ext_0;

	const struct const_Multiarray_R*const geom_coef = s_vol->geom_coef_p1;
	struct Multiarray_R* jacobian_vm = constructor_empty_Multiarray_R('C',3,(ptrdiff_t[]){n_vm,d,d}); // destructed
	for (ptrdiff_t row = 0; row < d; ++row)
		mm_NN1C_Operator_Multiarray_R(cv1_vg_vm.data[row],geom_coef,jacobian_vm,'d',2,NULL,&row);

	compute_cofactors_T((struct const_Multiarray_R*)jacobian_vm,(struct Multiarray_R*)s_vol->metrics_vm_p1);
	destructor_Multiarray_R(jacobian_vm);
}

static void compute_geometry_face_p1_T (struct Solver_Face_T* s_face)
{
	struct Face* face             = (struct Face*) s_face;
	struct Volume* vol            = face->neigh_info[0].volume;
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	struct const_Element* e            = (struct const_Element*) g_e;

	const int ind_lf = face->neigh_info[0].ind_lf;

	const struct Operator* vv0_vms_fgs = get_Multiarray_Operator(g_e->vv0_vms_fgs,(ptrdiff_t[]){ind_lf,0,0,1,1});

	const struct const_Multiarray_R* m_vm = s_vol->metrics_vm_p1;
	const struct const_Multiarray_R* metrics_p1 =
		constructor_mm_NN1_Operator_const_Multiarray_R(vv0_vms_fgs,m_vm,'C','d',m_vm->order,NULL); // d.

	compute_unit_normals_and_det_T(ind_lf,e->normals,metrics_p1,
		(struct Multiarray_R*)s_face->normals_p1,(struct Multiarray_R*)s_face->jacobian_det_p1);

	destructor_const_Multiarray_R(metrics_p1);
}

// Level 1 ********************************************************************************************************** //

/// \brief Correct the xyz coordinates of the internal face geometry nodes such that they are of degree 1.
static void correct_face_xyz_straight_T
	(struct Solver_Volume_T*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

static void compute_geom_coef_p1 (const struct Simulation*const sim, struct Solver_Volume_T*const s_vol)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;

	const struct Operator* vc0_vv_vgs = get_Multiarray_Operator(g_e->vc0_vv_vgs,(ptrdiff_t[]){0,0,1,1});

	const struct const_Multiarray_R*const xyz_ve = vol->xyz_ve;
	destructor_const_Multiarray_R(s_vol->geom_coef_p1);
	const_constructor_move_const_Multiarray_R(&s_vol->geom_coef_p1,
		constructor_mm_NN1_Operator_const_Multiarray_R(vc0_vv_vgs,xyz_ve,'C',op_format,xyz_ve->order,NULL)); // keep
}

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
	const struct const_Multiarray_R* xyz_s = constructor_xyz_s_ho_T('v',s_vol,sim); // destructed

	const struct const_Multiarray_R* xyz = constructor_xyz_blended_T('g',xyz_s,s_vol,sim); // destructed
	destructor_const_Multiarray_R(xyz_s);

	destructor_const_Multiarray_R(s_vol->geom_coef);
	const_constructor_move_const_Multiarray_R(&s_vol->geom_coef,constructor_geom_coef_ho_T(xyz,s_vol,sim)); // keep
	destructor_const_Multiarray_R(xyz);
}

static void compute_geom_coef_parametric_T (const struct Simulation*const sim, struct Solver_Volume_T*const s_vol)
{
	const struct const_Multiarray_R* xyz_s = constructor_xyz_s_ho_T('v',s_vol,sim); // destructed

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const struct const_Multiarray_R* xyz = test_case->constructor_xyz(0,xyz_s,s_vol,sim); // destructed
	destructor_const_Multiarray_R(xyz_s);

/// \todo add setup_bezier_mesh: Likely make this a separate function.

	destructor_const_Multiarray_R(s_vol->geom_coef);
	const_constructor_move_const_Multiarray_R(&s_vol->geom_coef,constructor_geom_coef_ho_T(xyz,s_vol,sim)); // keep
	destructor_const_Multiarray_R(xyz);

	if (is_internal_geom_straight()) {
		correct_face_xyz_straight_T(s_vol,sim);
		// Note: correct_internal_xyz_blended_T already called in correct_face_xyz_straight_T.
	} else {
		correct_internal_xyz_blended_T(s_vol,sim);
	}
}

// Level 2 ********************************************************************************************************** //

/** \brief Constructor for the 'f'ace 'geom'etry 'coef'ficients which are only accurate to degree p = 1.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_f_geom_coef_p1
	(const int side_index,                    ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face, ///< The \ref Solver_Face_T.
	 const struct Solver_Volume_T*const s_vol ///< The \ref Solver_Volume_T.
	);

/** \brief Return \ref constructor_const_Vector_i_inds_eq_1_const_Matrix_d  for the required
 *         \ref Geometry_Element::cc0_vgc_fgc operator.
 *  \return See brief. */
static const struct const_Vector_i* constructor_cc0_vgc_fgc_indices
	(const int side_index,                   ///< Side index of the volume whose geometry is being corrected.
	 const struct Solver_Face_T*const s_face ///< \ref Solver_Face_T.
	);

static void correct_face_xyz_straight_T (struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
{
	if (DIM == 1)
		return;

	const struct Volume*const vol = (struct Volume*) s_vol;
	for (int f  = 0; f  < NFMAX;    ++f) {
		const struct Face*const face = vol->faces[f][0];
		if (!face || face->boundary)
			continue;

		const int si = compute_side_index_face(face,vol);
		const struct Solver_Face_T*const s_face = (struct Solver_Face_T*) face;
		const struct const_Multiarray_R*const f_geom_coef_p1 = constructor_f_geom_coef_p1(si,s_face,s_vol); // dest.
		const struct const_Vector_i*const coef_inds = constructor_cc0_vgc_fgc_indices(si,s_face); // destructed

		update_rows_Multiarray_R((struct Multiarray_R*)s_vol->geom_coef,f_geom_coef_p1,coef_inds);

		destructor_const_Multiarray_R(f_geom_coef_p1);
		destructor_const_Vector_i(coef_inds);
	}

	correct_internal_xyz_blended_T(s_vol,sim);
}

// Level 3 ********************************************************************************************************** //

/** \brief Get the pointer to the appropriate \ref Geometry_Element::cv0_vgc_fis operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_vgc_fis
	(const int side_index,                   ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face ///< The \ref Solver_Face_T.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::vc0_fis_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__vc0_fis_fgc
	(const int side_index,                   ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face ///< The \ref Solver_Face_T.
	);

/** \brief Get the pointer to the appropriate \ref Geometry_Element::cc0_vgc_fgc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cc0_vgc_fgc
	(const int side_index,                   ///< The index of the side of the face under consideration.
	 const struct Solver_Face_T*const s_face ///< The current \ref Solver_Face_T.
	);

static const struct const_Multiarray_R* constructor_f_geom_coef_p1
	(const int side_index, const struct Solver_Face_T*const s_face, const struct Solver_Volume_T*const s_vol)
{
	const struct Operator* cv0_vgc_fis = get_operator__cv0_vgc_fis(side_index,s_face);
	const struct const_Multiarray_R*const geom_fis =
		constructor_mm_NN1_Operator_const_Multiarray_R(cv0_vgc_fis,s_vol->geom_coef,'C','d',2,NULL); // destructed

	const struct Operator* vc0_fis_fgc = get_operator__vc0_fis_fgc(side_index,s_face);
	const struct const_Multiarray_R*const geom_coef_fgc =
		constructor_mm_NN1_Operator_const_Multiarray_R(vc0_fis_fgc,geom_fis,'R','d',2,NULL); // returned
	destructor_const_Multiarray_R(geom_fis);

	return geom_coef_fgc;
}

static const struct const_Vector_i* constructor_cc0_vgc_fgc_indices
	(const int side_index, const struct Solver_Face_T*const s_face)
{
	const struct Operator*const cc0_vgc_fgc_op = get_operator__cc0_vgc_fgc(side_index,s_face);
	return constructor_const_Vector_i_inds_eq_1_const_Matrix_d(cc0_vgc_fgc_op->op_std);
}

// Level 4 ********************************************************************************************************** //

static const struct Operator* get_operator__cv0_vgc_fis (const int side_index, const struct Solver_Face_T*const s_face)
{
	const struct Face*const face             = (struct Face*) s_face;
	const struct Volume*const vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct Geometry_Element*const g_e  = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_lf = face->neigh_info[side_index].ind_lf,
	          p_v    = s_vol->p_ref;

	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->cv0_vgc_fis,(ptrdiff_t[]){ind_lf,0,0,1,p_v});
}

static const struct Operator* get_operator__vc0_fis_fgc (const int side_index, const struct Solver_Face_T*const s_face)
{
	const struct Face*const face             = (struct Face*) s_face;
	const struct Volume*const vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct Geometry_Element*const g_e  = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_e  = get_face_element_index(face),
	          p_v    = s_vol->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->vc0_fis_fgc,(ptrdiff_t[]){ind_e,ind_e,0,0,p_v,1});
}

static const struct Operator* get_operator__cc0_vgc_fgc (const int side_index, const struct Solver_Face_T*const s_face)
{
	const struct Face*const face             = (struct Face*) s_face;
	const struct Volume*const vol            = face->neigh_info[side_index].volume;
	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct Geometry_Element*const g_e  = &((struct Solver_Element*)vol->element)->g_e;

	const int ind_lf = face->neigh_info[side_index].ind_lf;
	const int p_v    = s_vol->p_ref;
	const int curved = ( (s_face->cub_type == 's') ? 0 : 1 );
	assert(curved);

	return get_Multiarray_Operator(g_e->cc0_vgc_fgc,(ptrdiff_t[]){ind_lf,0,0,p_v,p_v});
}
