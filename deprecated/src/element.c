// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "element.h"
#include "intrusive.h"

#include <string.h>

#include "constants_elements.h"
#include "Macros.h"
#include "multiarray.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for an individual \ref Element.
static struct Element* constructor_Element
	(const int elem_type ///< The element type.
	);

/// \brief Destructor for an individual \ref Element.
static void destructor_Element
	(struct Element* element ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct const_Intrusive_List* constructor_Element_List (const int d)
{
	struct Intrusive_List* Elements = constructor_empty_IL();

	push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(LINE));

	if (d >= 2) {
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(TRI));
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(QUAD));
	}

	if (d >= 3) {
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(TET));
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(HEX));
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(WEDGE));
		push_back_IL(Elements,(struct Intrusive_Link*) constructor_Element(PYR));
	}

	return (struct const_Intrusive_List*) Elements;
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
	int d,
	    n_f,
	    n_f_ve[NFMAX],
	    f_ve[NFMAX*NFVEMAX];
};

/** \brief Copy local (to each element type) element information to a container with larger scope.
 *	\returns Copy of local element info.
 */
static struct Elem_info copy_local_elem_info
	(const struct Elem_info*const src ///< The source element info.
	)
{
	struct Elem_info dest;

	dest.n_f = src->n_f;
	memcpy(dest.n_f_ve,src->n_f_ve,sizeof(src->n_f_ve));
	memcpy(dest.f_ve,  src->f_ve,  sizeof(src->f_ve));

	return dest;
}

static struct Element* constructor_Element
	(const int elem_type ///< The element type (e.g. LINE, TRI, ...)
	)
{
	struct Elem_info e_info;
	switch (elem_type) {
	case LINE: {
		const struct Elem_info e_info_l =
			{ .d      = 1,
			  .n_f    = 2,
			  .n_f_ve = {1, 1,},
			  .f_ve   = {0, 1,},
			};
		e_info = copy_local_elem_info(&e_info_l);
		break;
	} case TRI: {
		const struct Elem_info e_info_l =
			{ .d      = 2,
			  .n_f    = 3,
			  .n_f_ve = {2, 2, 2,},
			  .f_ve   = {1,2, 0,2, 0,1,},
			};
		e_info = copy_local_elem_info(&e_info_l);
		break;
	} case QUAD: {
		const struct Elem_info e_info_l =
			{ .d      = 2,
			  .n_f    = 4,
			  .n_f_ve = {2, 2, 2, 2,},
			  .f_ve   = {0,2, 1,3, 0,1, 2,3},
			};
		e_info = copy_local_elem_info(&e_info_l);
		break;
	} case TET: {
		EXIT_ADD_SUPPORT;
		break;
	} case HEX: {
		const struct Elem_info e_info_l =
			{ .d      = 3,
			  .n_f    = 6,
			  .n_f_ve = {4, 4, 4, 4, 4, 4,},
			  .f_ve   = {0,2,4,6, 1,3,5,7, 0,1,4,5, 2,3,6,7, 0,1,2,3, 4,5,6,7}
			};
		e_info = copy_local_elem_info(&e_info_l);
		break;
	} case WEDGE: {
		EXIT_ADD_SUPPORT;
		break;
	} case PYR: {
		EXIT_ADD_SUPPORT;
		break;
	} default: {
		EXIT_UNSUPPORTED;
		break;
	}}

	struct Multiarray_Vector_i* f_ve = constructor_copy_Multiarray_Vector_i(e_info.f_ve,e_info.n_f_ve,1,e_info.n_f);
//	print_Multiarray_Vector_i(f_ve);

	struct Element* element = malloc(sizeof *element); // returned

	const_cast_i(&element->type,elem_type);
	const_cast_i(&element->d,e_info.d);
	const_cast_i(&element->n_f,e_info.n_f);
	const_constructor_move_Multiarray_Vector_i(&element->f_ve,f_ve); // destructed


	return element;
}

static void destructor_Element (struct Element* element)
{
	destructor_Multiarray_Vector_i((struct Multiarray_Vector_i*)element->f_ve);
}
