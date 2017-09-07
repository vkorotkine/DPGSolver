// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "solver_volume.h"

#include "macros.h"
#include "constants_intrusive.h"

#include "simulation.h"
#include "geometry.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for an individual \ref Solver_Volume.
static struct Solver_Volume* constructor_Solver_Volume
	(struct Volume* volume ///< \ref Volume.
	);

/// \brief Destructor for an individual \ref Solver_Volume.
static void destructor_Solver_Volume
	(struct Solver_Volume* solver_volume ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_Solver_Volumes (struct Simulation*const sim)
{
	struct Intrusive_List* volumes        = sim->volumes;
	struct Intrusive_List* solver_volumes = constructor_empty_IL(IL_SOLVER_VOLUME);

	for (struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next)
		push_back_IL(solver_volumes,(struct Intrusive_Link*) constructor_Solver_Volume((struct Volume*) curr));

	set_up_geometry_solver(sim,solver_volumes);

EXIT_UNSUPPORTED;
	destructor_IL(volumes);
	return solver_volumes;
}

void destructor_Solver_Volumes (struct Intrusive_List* solver_volumes)
{
	for (const struct Intrusive_Link* curr = solver_volumes->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Solver_Volume((struct Solver_Volume*) curr);
		curr = next;
	}
	destructor_IL(solver_volumes);
}

void destructor_Solver_Volume (struct Solver_Volume* solver_volume)
{
// new members here.
	destructor_Volume((struct Volume*) solver_volume);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Solver_Volume* constructor_Solver_Volume (struct Volume* volume)
{
	struct Solver_Volume* solver_volume = calloc(1,sizeof *solver_volume); // returned

	solver_volume->volume = *volume; // shallow copy of the base.

// add new members here

	return solver_volume;
}
