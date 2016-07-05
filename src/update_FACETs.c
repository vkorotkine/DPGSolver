// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Update FACET information/operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	struct S_ELEMENT *ELEMENT;
}

void update_FACET_hp(void)
{
	// Initialize DB Parameters
	unsigned int Adapt = DB.Adapt;

	// Standard datatypes

	struct S_OPERATORS *OPS;
	struct S_FACET     *FACET;

	OPS = malloc(sizeof *OPS); // free

// Don't forget to free VOLUME->C_vC
	switch (Adapt) {
	default: // ADAPT_H, ADAPT_HP
		break;
	case ADAPT_P:
		// Change VIn/VOut so that FACET normals, etc are computed from VOLUME of higher order 
		break;
	}

	for (FACET = DB.FACET; FACET != NULL; FACET = FACET->next) {
//			init_ops(OPS,VOLUME,0);
	}
	free(OPS);
}
