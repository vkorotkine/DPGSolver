// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "Simulation.h"

#include <stdlib.h>

#include "Macros.h"


// Constructors/Destructors

struct Simulation* constructor_Simulation ()
{
	struct Simulation* sim = malloc(sizeof *sim); // returned;

	return sim;
}

void destructor_Simulation (struct Simulation* sim)
{
	FREE_NULL(sim);
}


// Setters/Getters

void set_Simulation_parameters
	(struct Simulation*const sim, unsigned int d, unsigned int n_var, unsigned int n_eq)
{
	*(unsigned int*)&sim->d     = d;
	*(unsigned int*)&sim->n_var = n_var;
	*(unsigned int*)&sim->n_eq  = n_eq;
}

void set_Simulation_element (struct Simulation*const sim, const struct Element*const e_head)
{
	*(const struct Element**)&sim->element_head = e_head;
}

void set_Simulation_volume (struct Simulation*const sim, struct Volume* v_head)
{
	sim->volume_head = v_head;
}

void set_Simulation_face (struct Simulation*const sim, struct Face* f_head)
{
	sim->face_head = f_head;
}
