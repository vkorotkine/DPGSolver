// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "geometry.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "constants_mesh.h"

#include "matrix.h"

#include "simulation.h"
#include "volume.h"
#include "intrusive.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to compute_geom_coef functions.
 *	\param sim    \ref Simulation.
 *	\param volume \ref Volume.
 */
typedef void (*compute_geom_coef_fptr)
	(const struct Simulation*const sim,
	 struct Volume*const volume
	);

/** \brief Set the appropriate function pointer for computing \ref Volume::geom_coef.
 *	\return See brief. */
static compute_geom_coef_fptr set_fptr_geom_coef
	(const int domain_type,   ///< \ref Simulation::domain_type.
	 const bool volume_curved ///< \ref Volume::curved.
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

void set_up_geometry_solver (struct Simulation* sim, struct Intrusive_List* volumes)
{
UNUSED(sim);
UNUSED(volumes);
//		compute_geom_metrics(sim,volume);
	EXIT_ADD_SUPPORT;
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

// Level 1 ********************************************************************************************************** //

static void compute_geom_coef_straight (const struct Simulation*const sim, struct Volume*const volume)
{
	destructor_Matrix_d((struct Matrix_d*)volume->geom_coef);

	if (strstr(sim->basis_geom,"lagrange") || strstr(sim->basis_geom,"bezier")) {
		const_constructor_copy_Matrix_d(&volume->geom_coef,volume->xyz_ve);
		if (volume->geom_coef->layout != 'C')
			transpose_Matrix_d((struct Matrix_d*)volume->geom_coef,true);
	} else if (strstr(sim->basis_geom,"nurbs")) {
		EXIT_ADD_SUPPORT;
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void compute_geom_coef_curved (const struct Simulation*const sim, struct Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
	EXIT_ADD_SUPPORT;
}

static void compute_geom_coef_parametric (const struct Simulation*const sim, struct Volume*const volume)
{
UNUSED(sim);
UNUSED(volume);
	EXIT_ADD_SUPPORT;
}
