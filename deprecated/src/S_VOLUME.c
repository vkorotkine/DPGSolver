// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "S_VOLUME.h"
#include "S_ELEMENT.h"

#include <stdlib.h>
#include <string.h>

/*
 *	Purpose:
 *		Provide functions related to the S_VOLUME struct.
 */


void set_element (struct S_VOLUME*const volume, struct S_ELEMENT*const element_head)
{
	/*
	 *	Comments:
	 *		'volume' must have been allocated dynamically for this to not be undefined behaviour.
	 *
	 *		The 'type' variable in volume should eventually be replaced with a pointer to its element (which holds the
	 *		type).
	 */

	*(struct S_ELEMENT**)&volume->element = get_element_by_type2(element_head,volume->type);
}
