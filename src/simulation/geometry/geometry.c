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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_intrusive.h"

#include "multiarray_operator.h"
#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "simulation.h"
#include "intrusive.h"
#include "solver_volume.h"
#include "geometry_element.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to compute_geom_coef functions.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Volume.
 */
typedef void (*compute_geom_coef_fptr)
	(const struct Simulation*const sim,
	 struct Volume*const volume
	);

/** \brief Set the appropriate function pointer for computing \ref Volume::geom_coef.
 *  \return See brief. */
static compute_geom_coef_fptr set_fptr_geom_coef
	(const int domain_type,   ///< \ref Simulation::domain_type.
	 const bool volume_curved ///< \ref Volume::curved.
	);

/** \brief Compute the geometry of the \ref Solver_Volume.
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
	(struct Simulation *sim,      ///< \ref Simulation.
	 struct Solver_Volume* volume ///< \ref Solver_Volume.
	);

// Interface functions ********************************************************************************************** //

void set_up_geometry (struct Simulation* sim, struct Intrusive_List* volumes)
{
	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Volume* volume = (struct Volume*) curr;

		compute_geom_coef_fptr compute_geom_coef = set_fptr_geom_coef(sim->domain_type,volume->curved);
		compute_geom_coef(sim,volume);
	}
}

void set_up_solver_geometry (struct Simulation* sim)
{
//	if ((sim->volumes->name != IL_SOLVER_VOLUME) || (sim->faces->name != IL_SOLVER_FACE))
	if ((sim->volumes->name != IL_SOLVER_VOLUME))
		EXIT_ERROR("Using incorrect volume (%d) and face (%d) lists.\n",sim->volumes->name,sim->faces->name);

	const struct const_Intrusive_List* geometry_elements = constructor_Geometry_Elements(sim); // destructed

	update_volumes_element(sim->volumes,geometry_elements);
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
//struct Volume* volume = (struct Volume*) curr;
//struct const_Geometry_Element* geometry_element = (struct const_Geometry_Element*) volume->element;
//printf("geom: %d %d %d\n",volume->index,((struct Element*)geometry_element)->type,volume->element->type);

		compute_geometry_volume(sim,(struct Solver_Volume*)curr);
	}
	update_volumes_element(sim->volumes,sim->elements);
EXIT_UNSUPPORTED;

	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
//		compute_geometry_face(sim,(struct Solver_Face*) curr);
	}

	destructor_const_IL(geometry_elements);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Compute \ref Volume::geom_coef for straight volumes.
static void compute_geom_coef_straight
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Volume*const volume         ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Volume::geom_coef for curved volumes using blending.
static void compute_geom_coef_curved
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Volume*const volume         ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Volume::geom_coef for curved volumes using the parametric mapping.
static void compute_geom_coef_parametric
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Volume*const volume         ///< Defined in \ref compute_geom_coef_fptr.
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

static void compute_geometry_volume (struct Simulation *sim, struct Solver_Volume* volume)
{
	struct Volume* base_volume = (struct Volume*) volume;
	struct const_Geometry_Element *element = (struct const_Geometry_Element*) base_volume->element;

	const int d = ((struct const_Element*)element)->d;

	const struct const_Multiarray_d*const geom_coef = base_volume->geom_coef;

	const int p = volume->p;

	struct Ops {
		const struct Multiarray_Operator* cv1_vg_vc;
	} ops =
		{ .cv1_vg_vc = constructor_default_Multiarray_Operator(),
		};

// Choose based on volume->curved
// Change to shorter name: "set_operator"
EXIT_ADD_SUPPORT; // Change to support using the Operator container
UNUSED(geom_coef);
UNUSED(d);
UNUSED(ops);
UNUSED(p);
/*
	set_const_Multiarray_Matrix_from_Multiarray_Matrix_d(ops.cv1_vg_vc,element->cv1_vgs_vcs,1,(ptrdiff_t[]){p,p,0});
//	set_const_Multiarray_Matrix_from_Multiarray_Matrix_d(ops.cv1_vg_vc,element->cv1_vgc_vcc,1,(ptrdiff_t[]){p,p,0});

	const ptrdiff_t n_vc = ops.cv1_vg_vc->data[0]->ext_0;

	struct Multiarray_d* jacobian_vc = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_vc,d,d});
//	struct Multiarray_d* jacobian_vm = constructor

	for (int row = 0; row < d; row++) {
		struct const_Vector_d geom_coef_V; // Note: not `const`
		set_const_Vector_from_Multiarray_d(&geom_coef_V,geom_coef,(ptrdiff_t[]){row});
		for (int col = 0; col < d; col++) {
			struct Vector_d jacobian_vc_V;
			set_Vector_from_Multiarray_d(&jacobian_vc_V,jacobian_vc,(ptrdiff_t[]){row,col});
			const struct const_Matrix_d* cv1_vg_vc = ops.cv1_vg_vc->data[col];

			mv_d('N',1.0,0.0,cv1_vg_vc,&geom_coef_V,&jacobian_vc_V);

//		mm_CTN_d(NvnI0,1,NvnG0,OPS->D_vG_vI[col],&XYZ[NvnG0*row],&J_vI[NvnI0*(d*row+col)]);
//		mm_CTN_d(NvnC0,1,NvnG0,OPS->D_vG_vC[col],&XYZ[NvnG0*row],&J_vC[NvnC0*(d*row+col)]);
		}
	}
*/

UNUSED(sim);
UNUSED(volume);
}

// Level 1 ********************************************************************************************************** //

static void compute_geom_coef_straight (const struct Simulation*const sim, struct Volume*const volume)
{
	destructor_const_Multiarray_d(volume->geom_coef);

	if (strstr(sim->basis_geom,"lagrange") || strstr(sim->basis_geom,"bezier")) {
		const_constructor_copy_Multiarray_d(&volume->geom_coef,volume->xyz_ve);
		if (volume->geom_coef->layout != 'C')
			transpose_Multiarray_d((struct Multiarray_d*)volume->geom_coef,true);
	} else if (strstr(sim->basis_geom,"nurbs")) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_ERROR("Unsupported sim->basis_geom: '%s'.",sim->basis_geom);
	}
}

static void compute_geom_coef_curved (const struct Simulation*const sim, struct Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
//	EXIT_ADD_SUPPORT;
}

static void compute_geom_coef_parametric (const struct Simulation*const sim, struct Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
	EXIT_ADD_SUPPORT;
}
