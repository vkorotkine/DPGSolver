// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_c.h"
#include "containers.h"
#include "containers_c.h"
#include "Simulation.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_VOLUME.h"

#include "compute_GradW_DG_c.h"

static struct Element_Solver_c* constructor_Elements_Solver_c (const struct Simulation*const simulation);
static struct Volume_Solver_c*  constructor_Volumes_Solver_c  (const struct Simulation*const simulation);
static void                     destructor_Elements_Solver_c  (struct Solver_c*const solver_c);
static void                     destructor_Volumes_Solver_c   (struct Solver_c*const solver_c);


struct Solver_c* constructor_Solver_c (const struct Simulation*const sim)
{
	struct Solver_c* solver_c = malloc(sizeof *solver_c); // returned

	// Parameters
	*(unsigned int*)&solver_c->d     = sim->d;
	*(unsigned int*)&solver_c->n_var = sim->n_var;

	// Elements
	*(struct Element_Solver_c**)&solver_c->element_head = constructor_Elements_Solver_c(sim); // destructed

	// Volumes
	solver_c->volume_head = constructor_Volumes_Solver_c(sim); // destructed

	return solver_c;
}

void destructor_Solver_c (struct Solver_c* solver_c)
{
	destructor_Elements_Solver_c(solver_c);
	destructor_Volumes_Solver_c(solver_c);
	FREE_NULL(solver_c);
}




// Static functions

// Element_Solver_c ************************************************************************************************* //

static struct Element_Solver_c* constructor_Element_Solver_c
	(const struct Element*const element)
{
}

static void destructor_Element_Solver_c (struct Element_Solver_c* element)
{
}

static struct Element_Solver_c* constructor_Elements_Solver_c
	(const struct Simulation*const simulation)
{
}

static void destructor_Elements_Solver_c
	(struct Solver_c*const solver_c)
{
}


// Volume_Solver_c ************************************************************************************************** //

static struct Volume_Solver_c* constructor_Volume_Solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME)
{
	struct Volume_Solver_c*const volume_s_c = calloc(1 , sizeof *volume_s_c); // returned

	const unsigned int P    = VOLUME->P,
	                   NvnG = VOLUME->NvnG,
	                   NvnS = VOLUME->NvnS,
	                   NvnI = { !VOLUME->curved ? VOLUME->element->NvnIs[P] : VOLUME->element->NvnIc[P] };

	const unsigned int d     = simulation->d,
	                   n_var = simulation->n_var;

	*(struct Multiarray_d**)&volume_s_c->xyz       = constructor_move_Multiarray_d_1_d('C',VOLUME->XYZ,2,NvnG,d);
	*(struct Multiarray_d**)&volume_s_c->det_jv_vi = constructor_move_Multiarray_d_1_d('C',VOLUME->detJV_vI,2,NvnI,1);
	*(struct Multiarray_d**)&volume_s_c->c_vi      = constructor_move_Multiarray_d_1_d('C',VOLUME->C_vI,3,NvnI,d,d);

	*(struct Multiarray_c**)&volume_s_c->w_hat   = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);
	*(struct Multiarray_c**)&volume_s_c->q_hat   = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d);
	*(struct Multiarray_c**)&volume_s_c->q_hat_v = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d);
	*(struct Multiarray_c**)&volume_s_c->rhs     = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);

	return volume_s_c;
}

static void destructor_Volume_Solver_c (struct Volume_Solver_c* volume)
{
	destructor_Multiarray_d_1((struct Multiarray_d*)volume->xyz);
	destructor_Multiarray_d_1((struct Multiarray_d*)volume->det_jv_vi);
	destructor_Multiarray_d_1((struct Multiarray_d*)volume->c_vi);

	destructor_Multiarray_c_1(volume->w_hat);
	destructor_Multiarray_c_1(volume->q_hat);
	destructor_Multiarray_c_1(volume->q_hat_v);
	destructor_Multiarray_c_1(volume->rhs);

	free(volume);
}

static struct Volume_Solver_c* constructor_Volumes_Solver_c (const struct Simulation*const simulation)
{
	struct S_VOLUME* VOLUME_head = simulation->volume_head;

	struct Volume_Solver_c* head = NULL,
	                      * prev = NULL;
	for (const struct S_VOLUME* VOLUME = VOLUME_head; VOLUME; VOLUME = VOLUME->next) {
		struct Volume_Solver_c* curr = constructor_Volume_Solver_c(simulation,VOLUME);

		if (prev)
			prev->next = curr;

		if (VOLUME == VOLUME_head)
			head = curr;

		prev = curr;
	}

	return head;
}

static void destructor_Volumes_Solver_c (struct Solver_c*const solver_c)
{
	struct Volume_Solver_c* next;
	for (struct Volume_Solver_c* volume = solver_c->volume_head; volume; ) {
		next = volume->next;
		destructor_Volume_Solver_c(volume);
		volume = next;
	}
}


void compute_GradW_c (const struct S_solver_info*const solver_info, const char stage)
{
	if (!(stage == 'A' || stage == 'F'))
		EXIT_UNSUPPORTED;

	switch (solver_info->method) {
	case METHOD_DG:
		if (stage == 'A')
			compute_GradW_DG_c(solver_info);
		else if (stage == 'F')
			free_GradW_DG_c(solver_info);

		break;
	case METHOD_HDG:
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}
