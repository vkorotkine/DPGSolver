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

#include "geometry.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "computational_elements.h"
#include "volume_solver.h"
#include "face_solver.h"

#include "const_cast.h"
#include "geometry_element.h"
#include "intrusive.h"
#include "operator.h"
#include "multiarray_operator.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

/** \brief Compute the geometry of the \ref Solver_Volume.
 *
 *  The following members are set:
 *  - Solver_Volume::metrics_vm;
 *  - Solver_Volume::metrics_vc;
 *  - Solver_Volume::jacobian_det_vc.
 *
 *  Following the analysis of Kopriva \cite Kopriva2006, the metric terms are computed using the curl-form such that
 *  the free-stream preservation property may be recovered. The consistent symmetric-conservative (CSC) metric of Abe
 *  and Haga (section 5.3) is used for the implementation \cite Abe2015. The steps are repeated below to clarify the
 *  equivalence of the prodecure adopted here with their procedure:
 *  - (step 0-1) As the metric contributions computed in step 0 are computed in a basis of sufficient order to
 *    represent them exactly and are subsequently interpolated to the consistent grid points (CGPs), the metric
 *    contributions are here computed directly at the CGPs.
 *  	- Our terminology for the GPs is R_vg ((R)eference coordinates of the (v)olume (g)eometry nodes).
 *  	- Our terminology for the CGPs is R_vm ((R)eference coordinates of the (v)olume (m)etric nodes).
 *  	- We allow for flexibility in the order of the R_vm nodes such that superparametric geometry can be used on
 *  	  curved domain boundaries; Abe and Haga use an isoparametric partial metric representation **before** the
 *  	  differentiation is applied, resulting in a subparametric metric representation (see eq. (43) \cite Abe2015).
 *  - (step 2) The computed metric terms are interpolated to the solution points (SPs).
 *  	- As the flux reconstruction scheme is collocated (solution interpolation and cubature nodes are coincident),
 *  	  the interpolation to the SPs is equivalent to the interpolation to the cubature nodes. Thus, interpolation
 *  	  to the R_vc ((R)eference coordinates of the (v)olume (c)ubature) is then performed in the implementation
 *  	  here.
 *
 *  \todo Investigate requirement of superparametric geometry on curved surfaces and add comments. Potentially ok by
 *        using over-integration in curved elements.
 *
 *  Given the 3D geometry Jacobian ordering of
 *
 *  \f{eqnarray*}{
 *  	J  = \{ &\{x_r,x_s,x_t\}, &\\
 *  	        &\{y_r,y_s,y_t\}, &\\
 *  	        &\{z_r,z_s,z_t\}  &\},
 *  \f}
 *
 *  using the nonconservative metric (NC) for clarity of exposition (section 5.1 \cite Abe2015), the ordering of the
 *  metric terms is:
 *
 *  \f{eqnarray*}{
 *  	m  = \{ &\{ +(y_s z_t - y_t z_s), -(y_r z_t - y_t z_r), +(y_r z_s - y_s z_r) \}, &\\
 *  	        &\{ -(x_s z_t - x_t z_s), +(x_r z_t - x_t z_r), -(x_r z_s - x_s z_r) \}, &\\
 *  	        &\{ +(x_s y_t - x_t y_s), -(x_r y_t - x_t y_r), +(x_r y_s - x_s y_r) \}, &\}.
 *  \f}
 */
static void compute_geometry_volume
	(struct Simulation* sim,      ///< \ref Simulation.
	 struct Solver_Volume* volume ///< \ref Solver_Volume.
	);

/** \brief Compute the geometry of the \ref Solver_Face.
 *
 *  The following members are set:
 *  - Solver_Face::xyz_fc;
 *  - Solver_Face::n_fc;
 *  - Solver_Face::jacobian_det_fc.
 */
static void compute_geometry_face
	(struct Simulation* sim,  ///< \ref Simulation.
	 struct Solver_Face* face ///< \ref Solver_Face.
	);

/// \brief Compute the face normal vectors at the nodes corresponding to the given face metrics.
static void compute_normals
	(const int ind_lf,                             ///< Defined for \ref compute_unit_normals_and_det.
	 const struct const_Multiarray_d* normals_ref, ///< Defined for \ref compute_unit_normals_and_det.
	 const struct const_Multiarray_d* metrics_f,   ///< Defined for \ref compute_unit_normals_and_det.
	 struct Multiarray_d* normals_f                ///< Defined for \ref compute_unit_normals_and_det.
	);

// Interface functions ********************************************************************************************** //

void set_up_solver_geometry (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);

	constructor_derived_Elements(sim,IL_ELEMENT_GEOMETRY);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		compute_geometry_volume(sim,(struct Solver_Volume*) curr);

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next)
		compute_geometry_face(sim,(struct Solver_Face*) curr);

	destructor_derived_Elements(sim,IL_ELEMENT);
}

void compute_unit_normals
	(const int ind_lf, const struct const_Multiarray_d* normals_ref, const struct const_Multiarray_d* metrics_f,
	 struct Multiarray_d* normals_f)
{
	compute_normals(ind_lf,normals_ref,metrics_f,normals_f);
	normalize_Multiarray_d(normals_f,"L2",false,NULL);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Function pointer to compute_geom_coef functions.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Volume.
 */
typedef void (*compute_geom_coef_fptr)
	(const struct Simulation*const sim,
	 struct Solver_Volume*const volume
	);

/** \brief Set the appropriate function pointer for computing \ref Solver_Volume::geom_coef.
 *  \return See brief. */
static compute_geom_coef_fptr set_fptr_geom_coef
	(const int domain_type,   ///< \ref Simulation::domain_type.
	 const bool volume_curved ///< \ref Volume::curved.
	);

/** \brief See return.
 *  \return The permutation required for conversion to the standard Jacobian ordering from the transposed ordering. */
static const ptrdiff_t* set_jacobian_permutation
	(const int d ///< The dimension
	);

/// \brief Compute the determinant of the geometry mapping Jacobian.
static void compute_detJV
	(struct const_Multiarray_d* jacobian, ///< Multiarray containing the Jacobian terms
	 struct Multiarray_d* jacobian_det    ///< Multiarray to be set to contain the determinant.
	);

/// \brief Compute the cofactor matrix entries of the geometry mapping Jacobian (referred to as the metrics).
static void compute_cofactors
	(struct const_Multiarray_d* jacobian, ///< Multiarray containing the Jacobian terms
	 struct Multiarray_d* metrics         ///< Multiarray to be set to contain the metric terms.
	);

/** \brief Compute the face unit normal vectors at the nodes corresponding to the given face metrics.
 *  The l^2 norm of the initially un-normalized normal vector at each of the nodes is stored in `jacobian_det`. This is
 *  done in accordance with the definition (see (eq. (B.6), \cite Zwanenburg2016)).
 */
static void compute_unit_normals_and_det
	(const int ind_lf,                             ///< \ref Face::Neigh_Info::ind_lf in \ref Face.
	 const struct const_Multiarray_d* normals_ref, ///< \ref Element::normals.
	 const struct const_Multiarray_d* metrics_f,   /**< \ref Solver_Volume::metrics_vm interpolated to the face
	                                                *   nodes. */
	 struct Multiarray_d* normals_f,               ///< \ref Multiarray_d\* in which to store the face normals.
	 struct Multiarray_d* jacobian_det_f           /**< \ref Multiarray_d\* in which to store the face jacobian
	                                                *   determinants. */
	);

static void compute_geometry_volume (struct Simulation *sim, struct Solver_Volume* volume)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) volume;
	const struct Geometry_Element* element = (struct Geometry_Element*) base_volume->element;

	compute_geom_coef_fptr compute_geom_coef = set_fptr_geom_coef(sim->domain_type,base_volume->curved);
	compute_geom_coef(sim,volume);

	const int d = ((struct const_Element*)element)->d;

	const struct const_Multiarray_d*const geom_coef = volume->geom_coef;

	const int p = volume->p_ref;

	struct Ops {
		const struct Multiarray_Operator* cv1_vg_vm;
		const struct Multiarray_Operator* cv1_vg_vc;
		const struct Operator* vv0_vm_vc;
	} ops = { .cv1_vg_vm = constructor_default_Multiarray_Operator(), // free (only)
	          .cv1_vg_vc = constructor_default_Multiarray_Operator(), // free (only)
	          .vv0_vm_vc = NULL, };

	const bool curved = base_volume->curved;
	const int p_g = ( curved ? p : 1 );
	set_MO_from_MO(ops.cv1_vg_vc,element->cv1_vg_vc[curved],1,(ptrdiff_t[]){0,0,p,p_g});
	set_MO_from_MO(ops.cv1_vg_vm,element->cv1_vg_vm[curved],1,(ptrdiff_t[]){0,0,p_g,p_g});
	ops.vv0_vm_vc = get_Multiarray_Operator(element->vv0_vm_vc[curved],(ptrdiff_t[]){0,0,p,p_g});

	const ptrdiff_t n_vm = ops.cv1_vg_vm->data[0]->op_std->ext_0,
	                n_vc = ops.cv1_vg_vc->data[0]->op_std->ext_0;

	struct Multiarray_d* jacobian_vm = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_vm,d,d}), // destructed
	                   * jacobian_vc = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_vc,d,d}); // destructed

	for (ptrdiff_t row = 0; row < d; ++row) {
		mm_NN1C_Operator_Multiarray_d(ops.cv1_vg_vm->data[row],geom_coef,jacobian_vm,op_format,2,NULL,&row);
		mm_NN1C_Operator_Multiarray_d(ops.cv1_vg_vc->data[row],geom_coef,jacobian_vc,op_format,2,NULL,&row);
	}

	free((void*)ops.cv1_vg_vm);
	free((void*)ops.cv1_vg_vc);

	const ptrdiff_t* perm = set_jacobian_permutation(d);
	permute_Multiarray_d(jacobian_vc,perm,jacobian_vc->layout);

	compute_detJV((struct const_Multiarray_d*)jacobian_vc,(struct Multiarray_d*)volume->jacobian_det_vc);
	compute_cofactors((struct const_Multiarray_d*)jacobian_vm,(struct Multiarray_d*)volume->metrics_vm);

	destructor_Multiarray_d(jacobian_vm);
	destructor_Multiarray_d(jacobian_vc);

	const struct const_Multiarray_d* met_vm = volume->metrics_vm;

	resize_Multiarray_d((struct Multiarray_d*)volume->metrics_vc,3,(ptrdiff_t[]){n_vc,d,d});
	mm_NN1C_Operator_Multiarray_d(
		ops.vv0_vm_vc,met_vm,(struct Multiarray_d*)volume->metrics_vc,op_format,met_vm->order,NULL,NULL);
}

static void compute_geometry_face (struct Simulation *sim, struct Solver_Face* face)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	assert((face->cub_type == 's') || (face->cub_type == 'c'));

	struct Face* base_face     = (struct Face*) face;
	struct Volume* base_volume = base_face->neigh_info[0].volume;
	struct Solver_Volume* volume = (struct Solver_Volume*) base_volume;

	const struct Geometry_Element* element = (struct Geometry_Element*) base_volume->element;
	struct const_Element* base_element     = (struct const_Element*) element;

	struct Ops {
		const struct Operator* cv0_vg_fc;
		const struct Operator* vv0_vm_fc;
	} ops = { .cv0_vg_fc = NULL,
	          .vv0_vm_fc = NULL, };

	const int ind_lf = base_face->neigh_info[0].ind_lf;
	const int p_v = volume->p_ref,
	          p_f = face->p_ref;

	const int curved_f = (face->cub_type == 's' ? 0 : 1);
	if (!base_volume->curved) {
		ops.cv0_vg_fc = get_Multiarray_Operator(element->cv0_vgs_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,1});
		ops.vv0_vm_fc = get_Multiarray_Operator(element->vv0_vms_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,1});
	} else {
		ops.cv0_vg_fc = get_Multiarray_Operator(element->cv0_vgc_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,p_v});
		ops.vv0_vm_fc = get_Multiarray_Operator(element->vv0_vmc_fc[curved_f],(ptrdiff_t[]){ind_lf,0,0,p_f,p_v});
	}

	const struct const_Multiarray_d* g_coef = volume->geom_coef;
	destructor_const_Multiarray_d(face->xyz_fc);
	const_constructor_move_const_Multiarray_d(&face->xyz_fc,
		constructor_mm_NN1_Operator_const_Multiarray_d(ops.cv0_vg_fc,g_coef,'C',op_format,g_coef->order,NULL)); // keep

	const struct const_Multiarray_d* m_vm = volume->metrics_vm;
	const struct const_Multiarray_d* metrics_fc =
		constructor_mm_NN1_Operator_const_Multiarray_d(ops.vv0_vm_fc,m_vm,'C',op_format,m_vm->order,NULL); // destructed

	compute_unit_normals_and_det(ind_lf,base_element->normals,metrics_fc,
		(struct Multiarray_d*)face->normals_fc,(struct Multiarray_d*)face->jacobian_det_fc);

	destructor_const_Multiarray_d(metrics_fc);
}

static void compute_normals
	(const int ind_lf, const struct const_Multiarray_d* normals_ref, const struct const_Multiarray_d* metrics_f,
	 struct Multiarray_d* normals_f)
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

	resize_Multiarray_d(normals_f,order_n,(ptrdiff_t[]){n_vals,d});
	set_to_value_Multiarray_d(normals_f,0.0);

	double* normals_d = normals_f->data;
	const double* normal_ref = get_row_const_Multiarray_d(ind_lf,normals_ref);
	const double* metrics_d  = metrics_f->data;
	for (ptrdiff_t dim_0 = 0; dim_0 < d; ++dim_0) {
	for (ptrdiff_t n = 0; n < n_vals; ++n) {
		for (ptrdiff_t dim_1 = 0; dim_1 < d; ++dim_1)
			*normals_d += normal_ref[dim_1]*metrics_d[n_vals*(dim_0+d*dim_1)+n];
		++normals_d;
	}}
	normals_f->layout = 'C';

	transpose_Multiarray_d(normals_f,true);
}

// Level 1 ********************************************************************************************************** //

/// \brief Compute \ref Solver_Volume::geom_coef for straight volumes.
static void compute_geom_coef_straight
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Solver_Volume::geom_coef for curved volumes using blending.
static void compute_geom_coef_curved
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Solver_Volume::geom_coef for curved volumes using the parametric mapping.
static void compute_geom_coef_parametric
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

static compute_geom_coef_fptr set_fptr_geom_coef (const int domain_type, const bool volume_curved)
{
	if (domain_type == DOM_STRAIGHT) {
		return compute_geom_coef_straight;
	} else if (domain_type == DOM_CURVED) {
		if (!volume_curved)
			return compute_geom_coef_straight;
		else
			return compute_geom_coef_curved;
	} else if (domain_type == DOM_PARAMETRIC) {
		return compute_geom_coef_parametric;
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

static void compute_detJV (struct const_Multiarray_d* jacobian, struct Multiarray_d* jacobian_det)
{
	const int order_j  = 3,
	          order_dj = 1;

	const ptrdiff_t* exts_j = jacobian->extents;

	assert(jacobian_det->order == order_dj);
	assert(jacobian->order == order_j);
	assert(exts_j[1] == exts_j[2]);

	const ptrdiff_t n_vals = exts_j[0],
	                d      = exts_j[1];

	resize_Multiarray_d(jacobian_det,order_dj,(ptrdiff_t[]){n_vals});
	double* j_det = jacobian_det->data;

	switch (d) {
	case 1: {
		const double* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})];
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			j_det[i] = x_r[i];
		break;
	} case 2: {
		const double* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
		            * x_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,1})],
		            * y_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,0})],
		            * y_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,1})];
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			j_det[i] = x_r[i]*y_s[i]-x_s[i]*y_r[i];
		break;
	} case 3: {
		const double* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
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

static void compute_cofactors (struct const_Multiarray_d* jacobian, struct Multiarray_d* metrics)
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

	resize_Multiarray_d(metrics,order_m,(ptrdiff_t[]){n_vals,d,d});
	switch (d) {
	case 1:
		for (ptrdiff_t i = 0; i < n_vals; ++i)
			metrics->data[i] = 1.0;
		break;
	case 2: {
		const double* x_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,0})],
		            * x_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){0,1})],
		            * y_r = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,0})],
		            * y_s = &jacobian->data[compute_index_sub_container(order_j,1,exts_j,(ptrdiff_t[]){1,1})];
		double* m_00 = &metrics->data[compute_index_sub_container(order_m,1,exts_m,(ptrdiff_t[]){0,0})],
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
		EXIT_ADD_SUPPORT; // curl-form. See comments in \ref compute_geometry_volume.
		break;
	} default:
		EXIT_ERROR("Unsupported: %td\n",d);
		break;
	}
}

static void compute_unit_normals_and_det
	(const int ind_lf, const struct const_Multiarray_d* normals_ref, const struct const_Multiarray_d* metrics_f,
	 struct Multiarray_d* normals_f, struct Multiarray_d* jacobian_det_f)
{
	compute_normals(ind_lf,normals_ref,metrics_f,normals_f);
	normalize_Multiarray_d(normals_f,"L2",true,jacobian_det_f);
}

// Level 2 ********************************************************************************************************** //

static void compute_geom_coef_straight (const struct Simulation*const sim, struct Solver_Volume*const volume)
{
	struct Volume* base_volume = (struct Volume*) volume;
	destructor_const_Multiarray_d(volume->geom_coef);

	if (strstr(sim->basis_geom,"lagrange") || strstr(sim->basis_geom,"bezier")) {
		const_constructor_copy_Multiarray_d(&volume->geom_coef,base_volume->xyz_ve);
		if (volume->geom_coef->layout != 'C')
			transpose_Multiarray_d((struct Multiarray_d*)volume->geom_coef,true);
	} else if (strstr(sim->basis_geom,"nurbs")) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_ERROR("Unsupported sim->basis_geom: '%s'.",sim->basis_geom);
	}
}

static void compute_geom_coef_curved (const struct Simulation*const sim, struct Solver_Volume*const volume)
{
	UNUSED(sim);
	struct Volume* base_volume = (struct Volume*) volume;
	const struct Geometry_Element* element = (struct Geometry_Element*) base_volume->element;

	const int p = volume->p_ref;

	const struct Operator* vc0_vg_vg = get_Multiarray_Operator(element->vc0_vgc_vgc,(ptrdiff_t[]){0,0,p,1});
// Potentially not the correct operator. vc0_vgs_vgc? vc0_vvs_vgc?

	const struct const_Multiarray_d* geom_coef =
		constructor_mm_NN1_Operator_const_Multiarray_d(vc0_vg_vg,base_volume->xyz_ve,'C','d',2,NULL); // keep
EXIT_ERROR("Add support after output to paraview is working.");

	destructor_const_Multiarray_d(volume->geom_coef);
	const_constructor_move_const_Multiarray_d(&volume->geom_coef,geom_coef);
}

static void compute_geom_coef_parametric (const struct Simulation*const sim, struct Solver_Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
	EXIT_ADD_SUPPORT;
}
