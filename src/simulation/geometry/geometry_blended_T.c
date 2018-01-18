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
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_tol.h"


#include "def_templates_geometry_blended.h"

#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref constructor_xyz_fptr_T used for the blended curved surface geometry corrections.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type,                     ///< The surface computational element type.
	 const char n_type,                      ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct const_Multiarray_R* xyz_i, ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Solver_Volume_T* s_vol,    ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Simulation* sim            ///< Defined for \ref constructor_xyz_fptr_T.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_R* constructor_xyz_blended_T
	(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const struct const_Multiarray_R*const xyz_e =
		( DIM == DMAX ? constructor_xyz_blended_ce('e',n_type,xyz_i,s_vol,sim) : xyz_i ); // dest. (if required)

	const struct const_Multiarray_R*const xyz = constructor_xyz_blended_ce('f',n_type,xyz_e,s_vol,sim); // returned
	if (DIM == DMAX)
		destructor_const_Multiarray_R(xyz_e);

	return xyz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for boundary computational element related data.
struct Boundary_Comp_Elem_Data {
	const struct const_Multiarray_R* xyz_ve; ///< \ref Volume::xyz_ve.

	int s_type; ///< \ref Element::s_type.

	const struct const_Vector_i* bc_boundaries;   ///< \ref Volume::bc_edges or \ref Volume::bc_faces.
	const struct const_Multiarray_Vector_i* b_ve; ///< \ref Element::e_ve or \ref Element::f_ve.

	/** See notation in \ref element_operators.h. The unknown parameter (X) may be set as:
	 *  - 'v'ertex (p2 vertices);
	 *  - 'g'eometry (standard blending).
	 */
	const struct Operator* vv0_vv_vX;

	const struct Multiarray_Operator* vv0_bv_vX; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.

	const struct Operator* vv0_vX_bX; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.
	const struct Operator* vv0_bX_vX; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.
	const struct Operator* vv0_vv_bX; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.
	const struct Operator* vv0_vv_bv; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.
	const struct Operator* vv0_bv_bX; ///< See \ref Boundary_Comp_Elem_Data::vv0_vv_vX.
};

/** \brief Constructor for a statically allocated \ref Boundary_Comp_Elem_Data containers.
 *  \return See brief. */
static struct Boundary_Comp_Elem_Data constructor_static_Boundary_Comp_Elem_Data
	(const char ce_type,                      ///< Computational element type. Options: 'e'dge, 'f'ace.
	 const char n_type,                       ///< Defined for \ref constructor_xyz_fptr_T.
	 const int p_geom,                        ///< The order of the geometry node basis.
	 const struct Solver_Volume_T*const s_vol ///< The current volume.
	);

/// \brief Destructor for a statically allocated \ref Boundary_Comp_Elem_Data containers.
static void destructor_static_Boundary_Comp_Elem_Data
	(struct Boundary_Comp_Elem_Data*const b_ce_d ///< Standard.
	);

/** \brief Compute the minimum polynomial degree to use as a base for the correction of the next approximation of the
 *         blended geometry.
 *  \return See brief.
 *
 *  Most commonly, the base degree is equal to the input volume polynomial degree. However, for certain blending
 *  functions, notably that of Lenoir (see \ref Simulation::geom_blending for reference), a sequential correction from
 *  the lowest (degree 2) to desired order is required and this function returns the lowest degree.
 */
static int compute_p_base_min
	(const struct Solver_Volume_T*const s_vol, ///< The current volume.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

/** \brief Set the operator members of \ref Boundary_Comp_Elem_Data which depend on the the current order, `p`, or
 *         boundary computational element index, `b`. */
static void set_Boundary_Comp_Elem_operators
	(struct Boundary_Comp_Elem_Data*const b_ce_d, ///< \ref Boundary_Comp_Elem_Data.
	 const struct Solver_Volume_T*const s_vol,    ///< The current volume.
	 const char ce_type,                          ///< The type of computational element.
	 const char n_type,                           ///< Defined for \ref constructor_xyz_fptr_T.
	 const int p,                                 ///< The current geometry correction order.
	 const int b                                  ///< The index of the boundary comp. element under consideration.
	);

/** \brief Constructor for the values of the blending to be used for each of the volume geometry nodes.
 *  \return See brief.
 *
 *  The blending functions for tensor-product and simplex element types are selected according to their definitions in
 *  \todo [add ref] Zwanenburg2017 (Discrete_Curvature).
 */
static const struct const_Vector_R* constructor_blend_values
	(const int ind_b,                                  ///< The index of the boundary under consideration.
	 const ptrdiff_t n_n,                              ///< The 'n'umber of geometry 'n'odes.
	 const struct Boundary_Comp_Elem_Data*const b_ce_d ///< \ref Boundary_Comp_Elem_Data.
	);

/** \brief Constructor for a \ref const_Matrix_T\* storing the difference between the corrected and current geometry
 *         nodal values on the boundary.
 *  \return See brief.
 *
 *  \note The difference for each of the volume geometry nodes **after** projection to the boundary is computed.
 */
static const struct const_Matrix_R* constructor_xyz_diff
	(const struct Boundary_Comp_Elem_Data*const b_ce_d, ///< \ref Boundary_Comp_Elem_Data.
	 const struct const_Matrix_R*const xyz_i,           ///< Input xyz coordinates.
	 const struct Solver_Volume_T*const s_vol           ///< The current \ref Solver_Volume_T.
	);

static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type, const char n_type, const struct const_Multiarray_R* xyz_i,
	 const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	assert(n_type == 'v' || n_type == 'g');
	const int p_geom = ( n_type == 'v' ? 2 : s_vol->p_ref ),
	          p_min = GSL_MIN(compute_p_base_min(s_vol,sim),p_geom);
	struct Boundary_Comp_Elem_Data b_ce_d =
		constructor_static_Boundary_Comp_Elem_Data(ce_type,n_type,p_geom,s_vol); // destructed

	struct Multiarray_R* xyz = constructor_copy_Multiarray_R((struct Multiarray_R*)xyz_i); // returned
	struct Matrix_R xyz_M = interpret_Multiarray_as_Matrix_R(xyz);

	const ptrdiff_t n_n = xyz_i->extents[0];
	for (int p = p_min; p <= p_geom; ++p) {
	for (int b = 0; b < b_ce_d.bc_boundaries->ext_0; ++b) {
		if (!is_bc_curved(b_ce_d.bc_boundaries->data[b]))
			continue;

		set_Boundary_Comp_Elem_operators(&b_ce_d,s_vol,ce_type,n_type,p,b);

		const struct const_Vector_R*const blend_values = constructor_blend_values(b,n_n,&b_ce_d); // destructed
		const struct const_Matrix_R*const xyz_diff =
			constructor_xyz_diff(&b_ce_d,(struct const_Matrix_R*)&xyz_M,s_vol); // destructed

		for (int d = 0; d < DIM; ++d) {
			Real*const data_xyz            = get_col_Matrix_R(d,&xyz_M);
			const Real*const data_xyz_diff = get_col_const_Matrix_R(d,xyz_diff);
			for (int n = 0; n < n_n; ++n)
				data_xyz[n] += blend_values->data[n]*data_xyz_diff[n];
		}

		destructor_const_Vector_d(blend_values);
		destructor_const_Matrix_R(xyz_diff);
	}}
	destructor_static_Boundary_Comp_Elem_Data(&b_ce_d);

	return (struct const_Multiarray_R*) xyz;
}

// Level 1 ********************************************************************************************************** //

static struct Boundary_Comp_Elem_Data constructor_static_Boundary_Comp_Elem_Data
	(const char ce_type, const char n_type, const int p_geom, const struct Solver_Volume_T*const s_vol)
{
	assert(DIM == DMAX || ce_type != 'e');
	assert(ce_type == 'e' || ce_type == 'f');

	struct Boundary_Comp_Elem_Data b_ce_d;

	struct Volume* vol = (struct Volume*) s_vol;
	const struct const_Element*const e = vol->element;

	b_ce_d.s_type = e->s_type;
	b_ce_d.xyz_ve = vol->xyz_ve;

	if (ce_type == 'f') {
		b_ce_d.bc_boundaries = vol->bc_faces;
		b_ce_d.b_ve          = e->f_ve;
	} else {
		b_ce_d.bc_boundaries = vol->bc_edges;
		b_ce_d.b_ve          = e->e_ve;
	}

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	if (n_type == 'g') {
		b_ce_d.vv0_vv_vX = get_Multiarray_Operator(g_e->vv0_vv_vg[1],(ptrdiff_t[]){0,0,p_geom,1});
		if (ce_type == 'f') {
			b_ce_d.vv0_bv_vX = alloc_MO_from_MO(g_e->vv0_fv_vgc,1,(ptrdiff_t[]){0,0,p_geom,1}); // free
		} else {
//			b_ce_d.vv0_bv_vX = set_MO_from_MO(g_e->vv0_ev_vg[1],1,(ptrdiff_t[]){0,0,p_geom,1});
			EXIT_ADD_SUPPORT; // Required for 3D.
		}
	} else {
//		b_ce_d.vv0_vv_vX = get_Multiarray_Operator(g_e->vv0_vv_vv,(ptrdiff_t[]){0,0,p_geom,1});
		EXIT_ADD_SUPPORT; // Required for h-adaptation.
	}

	return b_ce_d;
}

static void destructor_static_Boundary_Comp_Elem_Data (struct Boundary_Comp_Elem_Data*const b_ce_d)
{
	free_MO_from_MO(b_ce_d->vv0_bv_vX);
}

static int compute_p_base_min (const struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
{
	const struct Volume*const vol      = (struct Volume*) s_vol;
	const struct const_Element*const e = vol->element;

	int p_base_min = 0;

	const int s_type = e->s_type;
	switch (s_type) {
	case ST_TP:
/// \todo make `geom_blending` an array of `int`s such that switch statements can be used here.
		if (strcmp(sim->geom_blending[s_type],"gordon_hall") == 0)
			p_base_min = s_vol->p_ref;
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_blending[s_type]);
		break;
	case ST_SI:
		if (strcmp(sim->geom_blending[s_type],"szabo_babuska_gen") == 0 ||
		    strcmp(sim->geom_blending[s_type],"scott") == 0             ||
		    strcmp(sim->geom_blending[s_type],"lenoir_simple") == 0     ||
		    strcmp(sim->geom_blending[s_type],"nielson") == 0           )
			p_base_min = s_vol->p_ref;
		else if (strcmp(sim->geom_blending[s_type],"lenoir") == 0)
			p_base_min = 2;
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_blending[s_type]);
		break;
	case ST_PYR:
	case ST_WEDGE:
		EXIT_ADD_SUPPORT; ///< Ensure that all is working as expected.
	default:
		EXIT_ERROR("Unsupported: %d\n",s_type);
		break;
	}
	assert(p_base_min != 0);

	return p_base_min;
}

static void set_Boundary_Comp_Elem_operators
	(struct Boundary_Comp_Elem_Data*const b_ce_d, const struct Solver_Volume_T*const s_vol, const char ce_type,
	 const char n_type, const int p, const int b)
{
	assert(ce_type == 'e' || ce_type == 'f');

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	switch (n_type) {
	case 'g':
		if (ce_type == 'f') {
			b_ce_d->vv0_vX_bX = get_Multiarray_Operator(g_e->vv0_vgc_fgc,(ptrdiff_t[]){b,0,0,p,p});
			b_ce_d->vv0_bX_vX = get_Multiarray_Operator(g_e->vv0_fgc_vgc,(ptrdiff_t[]){b,0,0,p,p});
			b_ce_d->vv0_vv_bX = get_Multiarray_Operator(g_e->vv0_vv_fgc, (ptrdiff_t[]){b,0,0,p,1});
			b_ce_d->vv0_vv_bv = get_Multiarray_Operator(g_e->vv0_vv_fv,  (ptrdiff_t[]){b,0,0,1,1});
			b_ce_d->vv0_bv_bX = get_Multiarray_Operator(g_e->vv0_fv_fgc, (ptrdiff_t[]){0,0,0,0,p,1});
		} else {
			EXIT_ADD_SUPPORT;
		}
		break;
	case 'v':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_ERROR("Unsupported: %c",n_type);
		break;
	}
}

static const struct const_Vector_R* constructor_blend_values
	(const int ind_b, const ptrdiff_t n_n, const struct Boundary_Comp_Elem_Data*const b_ce_d)
{
	assert(DIM >= 2);
	struct Vector_d*const blend_values = constructor_empty_Vector_d(n_n); // returned
	double*const data_blend = blend_values->data;

	const struct Operator*const vv0_vv_vX = b_ce_d->vv0_vv_vX;
	const struct const_Vector_i*const b_ve_b = b_ce_d->b_ve->data[ind_b];

	switch (b_ce_d->s_type) {
	case ST_TP:
		for (int n = 0; n < n_n; ++n) {
			data_blend[n] = 0.0;
			const Real*const data_b_coords = get_row_const_Matrix_R(n,vv0_vv_vX->op_std);
			for (int ve = 0; ve < b_ve_b->ext_0; ++ve)
				data_blend[n] += data_b_coords[b_ve_b->data[ve]];
		}
		break;
	case ST_SI: {
		const struct Operator*const vv0_bv_vX = b_ce_d->vv0_bv_vX->data[ind_b];
		for (int n = 0; n < n_n; ++n) {
			const Real*const data_b_coords_num = get_row_const_Matrix_R(n,vv0_vv_vX->op_std),
			          *const data_b_coords_den = get_row_const_Matrix_R(n,vv0_bv_vX->op_std);

			Real blend_num = 1.0,
			     blend_den = 1.0;
			for (int ve = 0; ve < b_ve_b->ext_0; ++ve) {
				blend_num *= data_b_coords_num[b_ve_b->data[ve]];
				blend_den *= data_b_coords_den[ve];
			}
			data_blend[n] = ( (blend_num < EPS) ? 0.0 : blend_num/blend_den );
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",b_ce_d->s_type);
		break;
	}
	return (struct const_Vector_d*) blend_values;
}

static const struct const_Matrix_R* constructor_xyz_diff
	(const struct Boundary_Comp_Elem_Data*const b_ce_d, const struct const_Matrix_R*const xyz_i,
	 const struct Solver_Volume_T*const s_vol)
{
	const struct Volume*const vol = (struct Volume*) s_vol;
	struct Blended_Parametric_Data b_p_d =
		{ .xyz_ve = vol->xyz_ve,
		  .vv0_vv_bX = b_ce_d->vv0_vv_bX,
		  .vv0_vv_bv = b_ce_d->vv0_vv_bv,
		  .vv0_bv_bX = b_ce_d->vv0_bv_bX,
		  .normal = NULL,
		};

	const struct const_Matrix_R*const xyz_b = s_vol->constructor_xyz_surface(&b_p_d);

	struct Matrix_R*const xyz_surf_diff =
		constructor_mm_Matrix_R('N','N',-1.0,b_ce_d->vv0_vX_bX->op_std,xyz_i,'C'); // destructed
	add_in_place_Matrix_R(1.0,xyz_surf_diff,xyz_b);

	const struct const_Matrix_R*const xyz_diff = constructor_mm_const_Matrix_R(
		'N','N',1.0,b_ce_d->vv0_bX_vX->op_std,(struct const_Matrix_R*)xyz_surf_diff,'C'); // returned
	destructor_Matrix_R(xyz_surf_diff);

	return xyz_diff;
}
