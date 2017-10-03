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

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "const_cast.h"
#include "geometry_element.h"
#include "intrusive.h"
#include "operator.h"
#include "multiarray_operator.h"
#include "simulation.h"
#include "solver_volume.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to compute_geom_coef functions.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Volume.
 */
typedef void (*compute_geom_coef_fptr)
	(const struct Simulation*const sim,
	 struct Solver_Volume*const volume
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

void set_up_geometry (struct Simulation* sim, struct Intrusive_List* solver_volumes)
{
	EXIT_ERROR("Unused.\n"); /// \todo Delete this function. \deprecated
	for (struct Intrusive_Link* curr = solver_volumes->first; curr; curr = curr->next) {
		struct Volume*        volume        = (struct Volume*) curr;
		struct Solver_Volume* solver_volume = (struct Solver_Volume*) curr;

		compute_geom_coef_fptr compute_geom_coef = set_fptr_geom_coef(sim->domain_type,volume->curved);
		compute_geom_coef(sim,solver_volume);
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
		struct Volume* volume               = (struct Volume*) curr;
		struct Solver_Volume* solver_volume = (struct Solver_Volume*) curr;

		compute_geom_coef_fptr compute_geom_coef = set_fptr_geom_coef(sim->domain_type,volume->curved);
		compute_geom_coef(sim,solver_volume);

		compute_geometry_volume(sim,solver_volume);
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
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Volume::geom_coef for curved volumes using blending.
static void compute_geom_coef_curved
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Compute \ref Volume::geom_coef for curved volumes using the parametric mapping.
static void compute_geom_coef_parametric
	(const struct Simulation*const sim, ///< Defined in \ref compute_geom_coef_fptr.
	 struct Solver_Volume*const volume  ///< Defined in \ref compute_geom_coef_fptr.
	);

/// \brief Set the permutation required for conversion to the standard Jacobian ordering from the transposed ordering.
static const ptrdiff_t* set_jacobian_permutation
	(const int d ///< The dimension
	);

static compute_geom_coef_fptr set_fptr_geom_coef (const int domain_type, const bool volume_curved)
{
	if (domain_type == DOM_STRAIGHT) {
		return compute_geom_coef_straight;
	} else if (domain_type == DOM_CURVED) {
printf("vc: %d\n",volume_curved);
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
	struct const_Geometry_Element* element = (struct const_Geometry_Element*) base_volume->element;

	const int d = ((struct const_Element*)element)->d;

	const struct const_Multiarray_d*const geom_coef = volume->geom_coef;

	const int p = volume->p_ref;

	struct Ops {
		const struct Multiarray_Operator* cv1_vg_vc;
	} ops =
		{ .cv1_vg_vc = constructor_default_Multiarray_Operator(),
		};

	if (!base_volume->curved) {
		set_MO_from_MO(ops.cv1_vg_vc,element->cv1_vgs_vcs,1,(ptrdiff_t[]){0,0,1,1});
	} else {
		set_MO_from_MO(ops.cv1_vg_vc,element->cv1_vgc_vcc,1,(ptrdiff_t[]){0,0,p,p});
	}
print_Multiarray_Operator(ops.cv1_vg_vc);
print_const_Multiarray_d(geom_coef);

	const ptrdiff_t n_vc = ops.cv1_vg_vc->data[0]->op_std->ext_0;

	struct Multiarray_d* jacobian_vc = constructor_empty_Multiarray_d('C',3,(ptrdiff_t[]){n_vc,d,d});

	for (ptrdiff_t row = 0; row < d; ++row)
		mm_NN1C_Operator_Multiarray_d(ops.cv1_vg_vc->data[row],geom_coef,jacobian_vc,'d',2,NULL,&row);

	const ptrdiff_t* perm = set_jacobian_permutation(d);
print_Multiarray_d(jacobian_vc);
	permute_Multiarray_d(jacobian_vc,perm);
print_Multiarray_d(jacobian_vc);

// Choose based on volume->curved
// Change to shorter name: "set_operator"
EXIT_ADD_SUPPORT; // Change to support using the Operator container
UNUSED(geom_coef);
/*
//		mm_CTN_d(NvnI0,1,NvnG0,OPS->D_vG_vI[col],&XYZ[NvnG0*row],&J_vI[NvnI0*(d*row+col)]);
//		mm_CTN_d(NvnC0,1,NvnG0,OPS->D_vG_vC[col],&XYZ[NvnG0*row],&J_vC[NvnC0*(d*row+col)]);
		}
	}
*/

UNUSED(sim);
UNUSED(volume);
}

// Level 1 ********************************************************************************************************** //

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
	struct const_Geometry_Element* element = (struct const_Geometry_Element*) base_volume->element;

	const int p = volume->p_ref;

	const struct Operator* vc0_vg_vg = calloc(1,sizeof *vc0_vg_vg); // free

	set_O_from_MO(vc0_vg_vg,element->vc0_vgc_vgc,(ptrdiff_t[]){0,0,p,1});

	const struct const_Multiarray_d* geom_coef =
		constructor_mm_NN1_Operator_const_Multiarray_d(vc0_vg_vg,base_volume->xyz_ve,'C','d',2,NULL); // keep
EXIT_ERROR("Add support after output to paraview is working.");

	destructor_const_Multiarray_d(volume->geom_coef);
	const_constructor_move_const_Multiarray_d(&volume->geom_coef,geom_coef);

	free((void*)vc0_vg_vg);
}

static void compute_geom_coef_parametric (const struct Simulation*const sim, struct Solver_Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
	EXIT_ADD_SUPPORT;
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
