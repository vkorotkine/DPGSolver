// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_c.h"
#include "containers.h"
#include "containers_c.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_VOLUME.h"

#include "compute_GradW_DG_c.h"

/*	Purpose:
 *		Provide solver related functions for complex datatypes.
 */

static struct Volume_solver_c* constructor_Volume_solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME)
{
	struct Volume_solver_c*const volume_s_c = calloc(1 , sizeof *volume_s_c); // returned

	const unsigned int P    = VOLUME->P,
	                   NvnG = VOLUME->NvnG,
	                   NvnS = VOLUME->NvnS,
	                   NvnI = { !VOLUME->curved ? VOLUME->element->NvnIs[P] : VOLUME->element->NvnIc[P] };

	const unsigned int d     = simulation->d,
	                   n_var = simulation->n_var;

	*(struct Multiarray_d**)&volume_s_c->XYZ      = constructor_move_Multiarray_d_1_d('C',VOLUME->XYZ,2,NvnG,d);
	*(struct Multiarray_d**)&volume_s_c->detJV_vI = constructor_move_Multiarray_d_1_d('C',VOLUME->detJV_vI,2,NvnI,1);
	*(struct Multiarray_d**)&volume_s_c->C_vI     = constructor_move_Multiarray_d_1_d('C',VOLUME->C_vI,3,NvnI,d,d);

	*(struct Multiarray_c**)&volume_s_c->What  = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);   // keep
	*(struct Multiarray_c**)&volume_s_c->Qhat  = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d); // keep
	*(struct Multiarray_c**)&volume_s_c->QhatV = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d); // keep
	*(struct Multiarray_c**)&volume_s_c->RHS   = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);   // keep

	return volume_s_c;
}

static void destructor_Volume_solver_c (struct Volume_solver_c* volume)
{
	destructor_moved_Multiarray_d_1(volume->XYZ);
	destructor_moved_Multiarray_d_1(volume->detJV_vI);
	destructor_moved_Multiarray_d_1(volume->C_vI);

	destructor_Multiarray_c_1(volume->What);
	destructor_Multiarray_c_1(volume->Qhat);
	destructor_Multiarray_c_1(volume->QhatV);
	destructor_Multiarray_c_1(volume->RHS);

	free(volume);
}

static struct Volume_solver_c* constructor_Volumes_solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME_head)
{
	struct Volume_solver_c* head = NULL;
	                        prev = NULL;
	for (const struct S_VOLUME* VOLUME = VOLUME_head; VOLUME; VOLUME = VOLUME->next) {
		struct Volume_solver_c* curr = constructor_Volume_solver_c(simulation,VOLUME);

		if (prev)
			prev->next = curr;

		if (VOLUME == VOLUME_head)
			head = curr;

		prev = curr;
	}

	return head;
}

static void destructor_Volumes_solver_c (struct Volume_solver_c*const volume_head)
{
	struct Volume_solver_c* next;
	for (struct Volume_solver_c* volume = volume_head; volume; ) {
		next = volume->next;
		destructor_Volume_solver_c(volume);
		volume = next;
	}
}

struct Context_solver_c constructor_Context_solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME_head)
{

	struct Volume_solver_c*const volume_head = constructor_Volumes_solver_c(simulation,VOLUME_head);


	struct Context_solver_c context = { .d     = simulation->d,
	                                    .n_var = simulation->n_var,
	                                    .volume_head = volume_head,
	                                  };

	return context;
}

void destructor_Context_solver_c (struct Context_solver_c* context)
{
	destructor_Volumes_solver_c(context->volume_head);
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
