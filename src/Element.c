// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "Element.h"

#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "Intrusive.h"
#include "Multiarray.h"

struct Element* constructor_Element (const unsigned int elem_type);


struct Intrusive_List* constructor_Element_List (const unsigned int d)
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

	return Elements;
}

// Static functions ************************************************************************************************* //

/// \brief Constructs an \ref Element.
struct Element* constructor_Element
	(const unsigned int elem_type ///< The element type (e.g. LINE, TRI, ...)
	)
{
	struct Element* element = malloc(sizeof *element); // keep

	unsigned int Nf = 0,
	             n_f_ve[NFMAX]       = {0},
	             f_ve[NFMAX*NFVEMAX] = {0};

	switch (elem_type) {
	case LINE: {
		Nf = 2;
		const unsigned int n_f_ve_l[] = {1, 1,};
		const unsigned int f_ve_l[]   = {0, 1,};
		memcpy(n_f_ve,n_f_ve_l,sizeof(n_f_ve_l));
		memcpy(f_ve,f_ve_l,sizeof(f_ve_l));
		break;
	} case TRI: {
		Nf = 3;
		const unsigned int n_f_ve_l[] = {2, 2, 2,};
		const unsigned int f_ve_l[]   = {1,2, 0,2, 0,1,};
		memcpy(n_f_ve,n_f_ve_l,sizeof(n_f_ve_l));
		memcpy(f_ve,f_ve_l,sizeof(f_ve_l));
		break;
	} default: {
		EXIT_UNSUPPORTED;
		break;
	}}

	struct Multiarray_Vector_ui* corr_f_ve = constructor_empty_Multiarray_Vector_ui(1,Nf);
	set_Multiarray_Vector_ui_ui(corr_f_ve,f_ve,n_f_ve);
	print_Multiarray_Vector_ui(corr_f_ve);


	return element;
}
