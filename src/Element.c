// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "Element.h"
#include "Intrusive.h"

#include <string.h>

#include "constants_elements.h"
#include "Macros.h"
#include "Multiarray.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for an individual \ref Element.
static struct Element* constructor_Element
	(const unsigned int elem_type ///< The element type.
	);

/// \brief Destructor for an individual \ref Element.
static void destructor_Element
	(struct Element* element ///< Standard.
	);

// Interface functions ********************************************************************************************** //

struct const_Intrusive_List* constructor_Element_List (const unsigned int d)
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

struct const_Element* get_element_by_type (const struct const_Intrusive_List*const elements, const unsigned int type)
{
	for (const struct Intrusive_Link* curr = elements->first; curr; curr = curr->next) {
		struct const_Element* element = (struct const_Element*) curr;
		if (element->type == type)
			return element;
	}
	printf("Could not find the element of type: %d.\n",type);
	EXIT_UNSUPPORTED;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for local element-related information.
struct Elem_info {
	unsigned int d,
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
	(const unsigned int elem_type ///< The element type (e.g. LINE, TRI, ...)
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

	struct Multiarray_Vector_ui* f_ve = constructor_copy_Multiarray_Vector_ui(e_info.f_ve,e_info.n_f_ve,1,e_info.n_f);
//	print_Multiarray_Vector_ui(f_ve);

	struct Element* element = malloc(sizeof *element); // returned

	const_cast_ui(&element->type,elem_type);
	const_cast_ui(&element->d,e_info.d);
	const_cast_ui(&element->n_f,e_info.n_f);
	const_constructor_move_Multiarray_Vector_ui(&element->f_ve,f_ve); // destructed


	return element;
}

static void destructor_Element (struct Element* element)
{
	destructor_Multiarray_Vector_ui((struct Multiarray_Vector_ui*)element->f_ve);
}
