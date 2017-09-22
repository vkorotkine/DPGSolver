// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "S_ELEMENT.h"

#include <stdlib.h>

#include "Macros.h"

/*
 *	Purpose:
 *		Provide functions related to the S_ELEMENT struct.
 */

struct S_ELEMENT* get_element_by_type2 (struct S_ELEMENT*const element_head, const size_t type)
{
	for (struct S_ELEMENT* element = element_head; element; element = element->next) {
		if (element->type == type)
			return element;
	}
	EXIT_UNSUPPORTED;
}

const struct S_ELEMENT* get_const_element_by_type (const struct S_ELEMENT*const element_head, const size_t type)
{
	for (const struct S_ELEMENT* element = element_head; element; element = element->next) {
		if (element->type == type)
			return element;
	}
	EXIT_UNSUPPORTED;
}
