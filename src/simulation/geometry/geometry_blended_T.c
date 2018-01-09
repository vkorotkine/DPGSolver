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
#include "definitions_core.h"


#include "def_templates_geometry_blended.h"

#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //



/** \brief Version of \ref constructor_xyz_fptr_T used for the blended curved surface geometry corrections.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type,                     ///< The surface computational element type.
	 const struct const_Multiarray_R* xyz_i, ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Solver_Volume_T* s_vol,    ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Simulation* sim            ///< Defined for \ref constructor_xyz_fptr_T.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_R* constructor_xyz_blended_T
	(const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const struct const_Multiarray_R*const xyz_e =
		( DIM == DMAX ? constructor_xyz_blended_ce('e',xyz_i,s_vol,sim) : xyz_i ); // destructed (if required)

	const struct const_Multiarray_R*const xyz = constructor_xyz_blended_ce('f',xyz_e,s_vol,sim); // returned
	if (DIM == DMAX)
		destructor_const_Multiarray_R(xyz_e);

	return xyz;
}

const struct const_Multiarray_R* constructor_xyz_surface_cylinder_T
	(const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	UNUSED(xyz_i);
	UNUSED(s_vol);
	UNUSED(sim);

EXIT_UNSUPPORTED;
//	return (struct const_Multiarray_R*) xyz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for boundary computational element related data.
struct Boundary_Comp_Elem_Data {
	int n_b; ///< Number of boundary computational element entities.
};

/** \brief Constructor for a statically allocated \ref Boundary_Comp_Elem_Data containers.
 *  \return See brief. */
static struct Boundary_Comp_Elem_Data constructor_static_Boundary_Comp_Elem_Data
	(const char ce_type,           ///< Computational element type. Options: 'e'dge, 'f'ace.
	 const struct Volume*const vol ///< The current volume.
	);

/** \brief Compute the minimum polynomial degree to use as a base for the correction of the next approximation of the
 *         blended geometry.
 *  \return See brief.
 *
 *  Most commonly, the base degree is equal to the input volume polynomial degree. However, for certain blending
 *  functions, notably that of Lenoir (see \ref sim::geom_blending for reference), a sequential correction from the
 *  lowest (degree 2) to desired order is required and this function returns the lowest degree.
 */
static int compute_p_base_min
	(const struct Solver_Volume_T*const s_vol, ///< The current volume.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	struct Volume* vol = (struct Volume*) s_vol;
	const struct Boundary_Comp_Elem_Data b_ce_d = constructor_static_Boundary_Comp_Elem_Data(ce_type,vol);

	const ptrdiff_t n_n = xyz_i->extents[0];
	struct Multiarray_R* xyz = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
UNUSED(g_e);
	int p_geom = 0;
	if (!s_vol->computing_xyz_ve_p2) {
		p_geom = 2;
	} else {
		p_geom = s_vol->p_ref;
	}

	const int p_min = GSL_MIN(compute_p_base_min(s_vol,sim),p_geom);
	for (int p = p_min; p <= p_geom; ++p) {
	for (int b = 0; b < b_ce_d.n_b; ++b) {
	}}

EXIT_UNSUPPORTED;
	return (struct const_Multiarray_R*) xyz;
}

// Level 1 ********************************************************************************************************** //

static struct Boundary_Comp_Elem_Data constructor_static_Boundary_Comp_Elem_Data
	(const char ce_type, const struct Volume*const vol)
{
	assert(DIM == DMAX || ce_type != 'e');

	struct Boundary_Comp_Elem_Data b_ce_d;

	const struct const_Element*const e = vol->element;

	if (ce_type == 'e') {
		b_ce_d.n_b = e->n_e;
	} else if (ce_type == 'f') {
		b_ce_d.n_b = e->n_f;
	} else {
		EXIT_ERROR("Unsupported: %c\n",ce_type);
	}

	return b_ce_d;
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
