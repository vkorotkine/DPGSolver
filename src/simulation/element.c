// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "element.h"

#include <string.h>

#include "macros.h"
#include "constants_elements.h"
#include "constants_intrusive.h"

#include "multiarray.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for an individual \ref Element.
static struct Element* constructor_Element
	(const int elem_type ///< The element type (e.g. LINE, TRI, ...)
	);

// Interface functions ********************************************************************************************** //

struct const_Intrusive_List* constructor_Elements (const int d)
{
	struct Intrusive_List* elements = constructor_empty_IL(IL_ELEMENT);

	push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(LINE));

	if (d >= 2) {
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(TRI));
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(QUAD));
	}

	if (d >= 3) {
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(TET));
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(HEX));
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(WEDGE));
		push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(PYR));
	}

	return (struct const_Intrusive_List*) elements;
}

void destructor_Elements (struct Intrusive_List* elements)
{
	for (const struct Intrusive_Link* curr = elements->first; curr; ) {
		struct Intrusive_Link* next = curr->next;
		destructor_Element((struct Element*) curr);
		curr = next;
	}
	destructor_IL(elements);
}

void destructor_Element (struct Element* element)
{
	destructor_Multiarray_Vector_i((struct Multiarray_Vector_i*)element->f_ve);
}

void const_cast_const_Element (const struct const_Element*const* dest, const struct const_Element*const src)
{
	*(struct const_Element**) dest = (struct const_Element*) src;
}

struct const_Element* get_element_by_type (const struct const_Intrusive_List*const elements, const int type)
{
	for (const struct Intrusive_Link* curr = elements->first; curr; curr = curr->next) {
		struct const_Element* element = (struct const_Element*) curr;
		if (element->type == type)
			return element;
	}
	printf("Could not find the element of type: %d.\n",type);
	EXIT_UNSUPPORTED;
}

struct const_Element* get_element_by_face (const struct const_Element*const element, const int lf)
{
	int type_to_find = -1;
	switch (element->type) {
	case LINE:
		type_to_find = POINT;
		break;
	case TRI: case QUAD:
		type_to_find = LINE;
		break;
	case TET:
		type_to_find = TRI;
		break;
	case HEX:
		type_to_find = QUAD;
		break;
	case WEDGE:
		if (lf < 3)
			type_to_find = QUAD;
		else
			type_to_find = TRI;
		break;
	case PYR:
		if (lf < 4)
			type_to_find = TRI;
		else
			type_to_find = QUAD;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	for (const struct Intrusive_Link* curr = (const struct Intrusive_Link*) element; curr; curr = curr->prev) {
		struct const_Element* element_curr = (struct const_Element*) curr;
		if (element_curr->type == type_to_find)
			return element_curr;
	}
	EXIT_ERROR("Did not find the pointer to the face element");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for local element-related information.
struct Elem_info {
	int d,       ///< Defined in \ref Element.
	    n_ve,    ///< Defined in \ref Element.
	    n_f,     ///< Defined in \ref Element.
	    *n_f_ve, ///< The number of vertices on each face.
	    *f_ve;   ///< Defined in \ref Element.
};

static struct Element* constructor_Element (const int elem_type)
{
	// Note the use of the compound literals for the initialization of the local variables.
	struct Elem_info e_info;
	switch (elem_type) {
	case LINE:
		e_info.d    = 1;
		e_info.n_ve = 2;
		e_info.n_f  = 2;
		e_info.n_f_ve = (int[]) {1, 1,};
		e_info.f_ve   = (int[]) {0, 1,};
		break;
	case TRI:
		e_info.d    = 2;
		e_info.n_ve = 3;
		e_info.n_f  = 3;
		e_info.n_f_ve = (int[]) {2, 2, 2,};
		e_info.f_ve   = (int[]) {1,2, 0,2, 0,1,};
		break;
	case QUAD:
		e_info.d    = 2;
		e_info.n_ve = 4;
		e_info.n_f  = 4;
		e_info.n_f_ve = (int[]) {2, 2, 2, 2,};
		e_info.f_ve   = (int[]) {0,2, 1,3, 0,1, 2,3};
		break;
	case TET:
		EXIT_ADD_SUPPORT;
		break;
	case HEX:
		e_info.d    = 3;
		e_info.n_ve = 8;
		e_info.n_f  = 6;
		e_info.n_f_ve = (int[]) {4, 4, 4, 4, 4, 4,};
		e_info.f_ve   = (int[]) {0,2,4,6, 1,3,5,7, 0,1,4,5, 2,3,6,7, 0,1,2,3, 4,5,6,7};
		break;
	case WEDGE:
		EXIT_ADD_SUPPORT;
		break;
	case PYR:
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	const ptrdiff_t n_f = e_info.n_f;
	struct Multiarray_Vector_i* f_ve = constructor_copy_Multiarray_Vector_i_i(e_info.f_ve,e_info.n_f_ve,1,&n_f);

	struct Element* element = calloc(1,sizeof *element); // returned

	const_cast_i(&element->type,elem_type);
	const_cast_i(&element->d,e_info.d);
	const_cast_i(&element->n_ve,e_info.n_ve);
	const_cast_i(&element->n_f,e_info.n_f);
	const_constructor_move_Multiarray_Vector_i(&element->f_ve,f_ve); // destructed


	return element;
}
