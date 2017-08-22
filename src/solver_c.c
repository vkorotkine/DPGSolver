// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver_c.h"
#include "Multiarray.h"
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
	(const struct Simulation*const simulation, const struct S_ELEMENT*const element_base);
static void                     destructor_Element_Solver_c  (struct Element_Solver_c* element);

/// \brief Constructor for a doubly-linked list.
static struct Element_Solver_c* constructor_Elements_Solver_c
	(const struct Simulation*const simulation ///< Standard.
	)
{
	const struct S_ELEMENT*const element_head = (struct S_ELEMENT*) simulation->elements;

	struct Element_Solver_c* head = NULL,
	                       * prev = NULL;
	for (const struct S_ELEMENT* element = element_head; element; element = element->next) {
		struct Element_Solver_c* curr = constructor_Element_Solver_c(simulation,element);

		if (prev) {
			prev->next = curr;
			curr->prev = prev;
		}
		prev = curr;

		if (element == element_head)
			head = curr;
	}
	return head;
}

/// \brief Destructor for the list.
static void destructor_Elements_Solver_c
	(struct Solver_c*const solver_c ///< Standard.
	)
{
	for (struct Element_Solver_c* element = solver_c->element_head; element; ) {
		struct Element_Solver_c* next = element->next;
		destructor_Element_Solver_c(element);
		element = next;
	}
}
void set_up_operators_Solver_c (const struct Simulation*const simulation, const struct S_ELEMENT*const element_base);

/** \brief Constructor.
 *	Moves requires members of the base and constructs members which are necessary for the complex solver functions.
 */
static struct Element_Solver_c* constructor_Element_Solver_c
	(const struct Simulation*const simulation, ///< Standard.
	 const struct S_ELEMENT*const element_base   ///< The base \ref Element.
	)
{
	struct Element_Solver_c*const element = calloc(1 , sizeof *element); // returned

	//	Construct relevant operators
	set_up_operators_Solver_c(simulation,element_base);
	//	Move T
	//	Construct:
	//		- All cases: D_Strong
	//		- Not collocated: ChiS_vI, ChiS_fI,

	return element;
}

/// \brief Destructor.
static void destructor_Element_Solver_c
	(struct Element_Solver_c* element ///< Standard.
	)
{
UNUSED(element);
}

void set_up_operators_Solver_c
	(const struct Simulation*const simulation, ///< Standard.
	 const struct S_ELEMENT*const element_base   ///< The base \ref Element.
	)
{
	// E_vs_vc(s/c) (Previously ChiS_vI(s/c)) - Only using (s)traight in the comments below.
	// Needs Er_vs_vcs*T_vs_vs
	// Needs cub data, basis evaluation.
UNUSED(simulation); // Need simulation->collocated
UNUSED(element_base);
}


// Volume_Solver_c ************************************************************************************************** //
static struct Volume_Solver_c* constructor_Volume_Solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const volume_base);
static void destructor_Volume_Solver_c (struct Volume_Solver_c* volume);

/// \brief Constructor for a doubly-linked list.
static struct Volume_Solver_c* constructor_Volumes_Solver_c
	(const struct Simulation*const simulation ///< Standard.
	)
{
	struct S_VOLUME* volume_base_head = (struct S_VOLUME*) simulation->volumes;

	struct Volume_Solver_c* head = NULL,
	                      * prev = NULL;
	for (const struct S_VOLUME* volume_base = volume_base_head; volume_base; volume_base = volume_base->next) {
		struct Volume_Solver_c* curr = constructor_Volume_Solver_c(simulation,volume_base);

		if (prev) {
			curr->prev = prev;
			curr->prev->next = curr;
		}
		prev = curr;

		if (volume_base == volume_base_head)
			head = curr;
	}
	return head;
}

/// \brief Destructor for the list.
static void destructor_Volumes_Solver_c
	(struct Solver_c*const solver_c ///< Standard.
	)
{
	for (struct Volume_Solver_c* volume = solver_c->volume_head; volume; ) {
		struct Volume_Solver_c* next = volume->next;
		destructor_Volume_Solver_c(volume);
		volume = next;
	}
}

/** \brief Constructor.
 *	Moves requires members of the base and constructs members which are necessary for the complex solver functions.
 */
static struct Volume_Solver_c* constructor_Volume_Solver_c
	(const struct Simulation*const simulation, ///< Standard.
	 const struct S_VOLUME*const volume_base     ///< The base \ref Volume.
	)
{
	struct Volume_Solver_c*const volume = calloc(1 , sizeof *volume); // returned

	const unsigned int P    = volume_base->P,
	                   NvnG = volume_base->NvnG,
	                   NvnS = volume_base->NvnS,
	                   NvnI = (!volume_base->curved ? volume_base->element->NvnIs[P] : volume_base->element->NvnIc[P]);

	const unsigned int d     = simulation->d,
	                   n_var = simulation->n_var;

	*(struct Multiarray_d**)&volume->xyz       = constructor_move_Multiarray_d_d('C',volume_base->XYZ,2,NvnG,d);
	*(struct Multiarray_d**)&volume->det_jv_vi = constructor_move_Multiarray_d_d('C',volume_base->detJV_vI,2,NvnI,1);
	*(struct Multiarray_d**)&volume->c_vi      = constructor_move_Multiarray_d_d('C',volume_base->C_vI,3,NvnI,d,d);
	// if (!collocated) op_M_inv

	// Will likely require modification for HDG.
	*(struct Multiarray_c**)&volume->w_hat   = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);
	*(struct Multiarray_c**)&volume->q_hat   = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d);
	*(struct Multiarray_c**)&volume->q_hat_v = constructor_empty_Multiarray_c_1('C',3,NvnS,n_var,d);
	*(struct Multiarray_c**)&volume->rhs     = constructor_empty_Multiarray_c_1('C',2,NvnS,n_var);

	return volume;
}

/// \brief Destructor.
static void destructor_Volume_Solver_c
	(struct Volume_Solver_c* volume ///< Standard.
	)
{
	destructor_Multiarray_d((struct Multiarray_d*)volume->xyz);
	destructor_Multiarray_d((struct Multiarray_d*)volume->det_jv_vi);
	destructor_Multiarray_d((struct Multiarray_d*)volume->c_vi);

	destructor_Multiarray_c_1(volume->w_hat);
	destructor_Multiarray_c_1(volume->q_hat);
	destructor_Multiarray_c_1(volume->q_hat_v);
	destructor_Multiarray_c_1(volume->rhs);

	free(volume);
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
