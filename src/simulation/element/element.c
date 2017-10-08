/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/// \file

#include "element.h"

#include <assert.h>
#include <string.h>
#include <math.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_intrusive.h"
#include "definitions_math.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for an individual \ref Element.
 *  \return Standard. */
static struct Element* constructor_Element
	(const int elem_type ///< The element type (e.g. LINE, TRI, ...)
	);

/// \brief Set \ref Element::present to `true` for element of the input type and its associated elements.
void set_element_present
	(const int e_type,                           ///< \ref Element::type.
	 const struct const_Intrusive_List* elements ///< \ref Simulation::elements.
	);

/** \brief Mutable version of \ref get_element_by_type.
 *  \return See brief. */
static struct Element* get_mutable_element_by_type
	(const struct Intrusive_List*const elements, ///< Defined for \ref get_element_by_type.
	 const int type                              ///< Defined for \ref get_element_by_type.
	);

// Interface functions ********************************************************************************************** //

struct const_Intrusive_List* constructor_Elements (const int d)
{
	assert(sizeof(struct Element) == sizeof(struct const_Element));

	struct Intrusive_List* elements = constructor_empty_IL(IL_ELEMENT,NULL);

	push_back_IL(elements,(struct Intrusive_Link*) constructor_Element(POINT));
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

	set_tp_sub_elements(elements);

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

void destructor_const_Elements (const struct const_Intrusive_List* elements)
{
	destructor_Elements((struct Intrusive_List*)elements);
}

void destructor_Element (struct Element* element)
{
	destructor_Multiarray_Vector_i(element->f_ve);
	destructor_Multiarray_d(element->normals);
}

void const_cast_const_Element (const struct const_Element*const* dest, const struct const_Element*const src)
{
	*(struct const_Element**) dest = (struct const_Element*) src;
}

const struct const_Element* get_element_by_type (const struct const_Intrusive_List*const elements, const int type)
{
	for (const struct const_Intrusive_Link* curr = elements->first; curr; curr = curr->next) {
		struct const_Element* element = (struct const_Element*) curr;
		if (element->type == type)
			return element;
	}
	EXIT_ERROR("Could not find the element of type: %d.\n",type);
}

const struct const_Element* get_element_by_face (const struct const_Element*const element, const int lf)
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

bool wedges_present (const struct const_Intrusive_List*const elements)
{
	for (const struct const_Intrusive_Link* curr = elements->first; curr; curr = curr->next) {
		struct const_Element* element = (struct const_Element*) curr;
		if (element->type == WEDGE)
			return true;
	}
	return false;
}

int compute_elem_type_sub_ce (const int e_type, const char ce, const int ind_ce)
{
	switch (ce) {
	case 'v':
		switch (e_type) {
		case LINE: // fallthrough
		case TRI:  // fallthrough
		case QUAD: // fallthrough
		case TET:  // fallthrough
		case HEX:  // fallthrough
		case WEDGE:
			return e_type;
		case PYR:
			switch (ind_ce) {
			case 0: // standard pyr
			case 1: case 2:  case 3: case 4: // 4 base pyr sub-elements.
			case 9: case 10: // 2 top pyr sub-elements.
				return PYR;
				break;
			case 5: case 6: case 7: case 8: // 4 tet sub-elements.
				return TET;
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_ce);
				break;
			}
		default:
			EXIT_ERROR("Unsupported: %d\n",e_type);
			break;
		}
		break;
	case 'f':
		switch (e_type) {
		case LINE:
			return POINT;
			break;
		case TRI: // fallthrough
		case QUAD:
			return LINE;
			break;
		case TET:
			return TRI;
			break;
		case HEX:
			return QUAD;
			break;
		case WEDGE:
/// \todo Replace with definitions from definitions_h_ref.h
			switch (ind_ce) {
			case 0:  case 1:  case 2:           // 3 quad faces of wedge reference element.
			case 5:  case 6:  case 7:  case 8:  // 4 quad refined face 0 sub-faces.
			case 9:  case 10: case 11: case 12: // 4 quad refined face 1 sub-faces.
			case 13: case 14: case 15: case 16: // 4 quad refined face 2 sub-faces.
				return QUAD;
				break;
			case 3: case 4:                     // 2 tri faces of wedge reference element.
			case 17: case 18: case 19: case 20: // 4 tri refined face 3 sub-faces.
			case 21: case 22: case 23: case 24: // 4 tri refined face 4 sub-faces.
				return TRI;
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_ce);
				break;
			}
			break;
		case PYR:
			switch (ind_ce) {
			case 0:  case 1:  case 2:  case 3:  // 4 tri faces of pyr reference element.
			case 5:  case 6:  case 7:  case 8:  // 4 tri refined face 0 sub-faces.
			case 9:  case 10: case 11: case 12: // 4 tri refined face 1 sub-faces.
			case 13: case 14: case 15: case 16: // 4 tri refined face 2 sub-faces.
			case 17: case 18: case 19: case 20: // 4 tri refined face 3 sub-faces.
				return TRI;
				break;
			case 4:                             // 1 quad face of pyr reference element.
			case 21: case 22: case 23: case 24: // 4 quad refined face 4 sub-faces.
				return QUAD;
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_ce);
				break;
			}
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",e_type);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",ce);
		break;
	}
}

int compute_super_from_elem_type (const int e_type)
{
	switch (e_type) {
	case POINT: case LINE: case QUAD: case HEX:
		return ST_TP;
		break;
	case TRI: case TET:
		return ST_SI;
		break;
	case PYR:
		return ST_PYR;
		break;
	case WEDGE:
		return ST_WEDGE;
		break;
	default:
		EXIT_ERROR("Unsupported (%d)\n",e_type);
		break;
	}
}

int compute_elem_from_super_type (const int s_type, const int d)
{
	switch (s_type) {
	case ST_TP:
		switch (d) {
			case 0: return POINT; break;
			case 1: return LINE;  break;
			case 2: return QUAD;  break;
			case 3: return HEX;   break;
			default:
				EXIT_ERROR("Unsupported: %d\n",d);
				break;
		}
		break;
	case ST_SI:
		if      (d == 2) return TRI;
		else if (d == 3) return TET;
		else             EXIT_ERROR("Unsupported: %d\n",d);
		break;
	case ST_PYR:
		assert(d == 3);
		return PYR;
		break;
	case ST_WEDGE:
		assert(d == 3);
		return WEDGE;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",s_type);
		break;
	}
}

void set_elements_present (const struct const_Intrusive_List* elements, const struct const_Vector_i*const elem_types)
{
	int e_type = -1;

	const ptrdiff_t n_elem = elem_types->ext_0;
	for (ptrdiff_t i = 0; i < n_elem; ++i) {
		if (e_type != elem_types->data[i]) {
			e_type = elem_types->data[i];
			set_element_present(e_type,elements);
		}
	}
}

void set_tp_sub_elements (struct Intrusive_List* elements)
{
	for (const struct Intrusive_Link* curr = elements->first; curr; curr = curr->next) {
		struct Element* element = (struct Element*) curr;
		switch (element->type) {
		case POINT: // fallthrough
		case LINE:  // fallthrough
		case TRI:   // fallthrough
		case TET:   // fallthrough
		case PYR:
			element->sub_element[0] = NULL;
			element->sub_element[1] = NULL;
			break;
		case QUAD: // fallthrough
		case HEX:
			element->sub_element[0] = get_mutable_element_by_type(elements,LINE);
			element->sub_element[1] = get_mutable_element_by_type(elements,LINE);
			break;
		case WEDGE:
			element->sub_element[0] = get_mutable_element_by_type(elements,TRI);
			element->sub_element[1] = get_mutable_element_by_type(elements,LINE);
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",element->type);
			break;
		}
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for local element-related information.
struct Elem_info {
	int s_type,  ///< Defined in \ref Element.
	    d,       ///< Defined in \ref Element.
	    n_ve,    ///< Defined in \ref Element.
	    n_e,     ///< Defined in \ref Element.
	    n_f,     ///< Defined in \ref Element.
	    *n_f_ve, ///< The number of vertices on each face.
	    *f_ve;   ///< Defined in \ref Element.

	int n_ref_max_v, ///< Defined in \ref Element.
	    n_ref_max_f; ///< Defined in \ref Element.
};

/** \brief Constructor for \ref Element::normals.
 *  \return See brief.
 *
 *  The normals are computed using (eq. (B.6), \cite Zwanenburg2016).
 */
struct Multiarray_d* constructor_reference_normals
	(const int e_type,              ///< The element type.
	 const struct Elem_info* e_info ///< \ref Elem_info.
	);

static struct Element* constructor_Element (const int elem_type)
{
	struct Elem_info e_info;
	switch (elem_type) {
	case POINT:
		e_info.s_type = ST_TP;
		e_info.d      = 0;
		e_info.n_ve   = 1;
		e_info.n_e    = 0;
		e_info.n_f    = 0;
		e_info.n_f_ve = NULL;
		e_info.f_ve   = NULL;

		e_info.n_ref_max_v = 1;
		e_info.n_ref_max_f = 0;
		break;
	case LINE:
		e_info.s_type = ST_TP;
		e_info.d      = 1;
		e_info.n_ve   = 2;
		e_info.n_e    = 2;
		e_info.n_f    = 2;
		e_info.n_f_ve = (int[]) {1, 1,};
		e_info.f_ve   = (int[]) {0, 1,};

		e_info.n_ref_max_v = 3;
		e_info.n_ref_max_f = 1;
		break;
	case TRI:
		e_info.s_type = ST_SI;
		e_info.d      = 2;
		e_info.n_ve   = 3;
		e_info.n_e    = 3;
		e_info.n_f    = 3;
		e_info.n_f_ve = (int[]) {2, 2, 2,};
		e_info.f_ve   = (int[]) {1,2, 0,2, 0,1,};

		e_info.n_ref_max_v = 5;
		e_info.n_ref_max_f = 3;
		break;
	case QUAD:
		e_info.s_type = ST_TP;
		e_info.d      = 2;
		e_info.n_ve   = 4;
		e_info.n_e    = 4;
		e_info.n_f    = 4;
		e_info.n_f_ve = (int[]) {2, 2, 2, 2,};
		e_info.f_ve   = (int[]) {0,2, 1,3, 0,1, 2,3};

		e_info.n_ref_max_v = 5;
		e_info.n_ref_max_f = 3;
		break;
	case TET:
		EXIT_ADD_SUPPORT;
		break;
	case HEX:
		e_info.s_type = ST_TP;
		e_info.d      = 3;
		e_info.n_ve   = 8;
		e_info.n_e    = 12;
		e_info.n_f    = 6;
		e_info.n_f_ve = (int[]) {4, 4, 4, 4, 4, 4,};
		e_info.f_ve   = (int[]) {0,2,4,6, 1,3,5,7, 0,1,4,5, 2,3,6,7, 0,1,2,3, 4,5,6,7};

		e_info.n_ref_max_v = 9;
		e_info.n_ref_max_f = 5;
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

	struct Element* element = calloc(1,sizeof *element); // returned

	element->present = false;
	element->type = elem_type;
	element->s_type = e_info.s_type;
	element->d = e_info.d;
	element->n_ve = e_info.n_ve;
	element->n_e = e_info.n_e;
	element->n_f = e_info.n_f;
	element->n_ref_max_v = e_info.n_ref_max_v;
	element->n_ref_max_f = e_info.n_ref_max_f;
	element->f_ve = constructor_copy_Multiarray_Vector_i_i(e_info.f_ve,e_info.n_f_ve,1,&n_f); // destructed

	element->normals = constructor_reference_normals(elem_type,&e_info); // destructed

	element->derived = NULL;

	return element;
}

void set_element_present (const int e_type, const struct const_Intrusive_List* elements)
{
	// Enable if needed.
	if (e_type == POINT)
		return;

	struct Element* element = (struct Element*) get_element_by_type(elements,e_type);
	element->present = true;

	switch (e_type) {
	case POINT:
		// Do nothing.
		break;
	case LINE:
		set_element_present(POINT,elements); // Enable if needed.
		break;
	case TRI: // fallthrough
	case QUAD:
		set_element_present(LINE,elements);
		break;
	case TET:
		set_element_present(TRI,elements);
		break;
	case HEX:
		set_element_present(QUAD,elements);
		break;
	case WEDGE:
		set_element_present(LINE,elements);
		set_element_present(TRI,elements);
//		set_element_present(QUAD,elements); // Should not be necessary.
		break;
	case PYR:
		set_element_present(TRI,elements);
		set_element_present(QUAD,elements);
		set_element_present(TET,elements); // For h-refinement.
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
}

static struct Element* get_mutable_element_by_type (const struct Intrusive_List*const elements, const int type)
{
	return (struct Element*) get_element_by_type((const struct const_Intrusive_List*)elements,type);
}

// Level 1 ********************************************************************************************************** //

///\{ \name Constants related to reference normals.
#define THETA_E_TET 1.230959417340774682134929178 // atan(2.0*sqrt(2.0))
#define THETA_E_PYR 0.955316618124509278163857103 // atan(sqrt(2.0))
///\}

struct Multiarray_d* constructor_reference_normals (const int e_type, const struct Elem_info* e_info)
{
	const double* theta_eta  = NULL,
	            * theta_zeta = NULL;

	switch (e_type) {
	case POINT:
		theta_eta  = NULL;
		theta_zeta = NULL;
		break;
	case LINE: {
		static const double t_e[] = { 0.0, 0.0, };
		static const double t_z[] = { PI, 0.0, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case TRI: {
		static const double t_e[] = { 0.0, 0.0, 0.0, };
		static const double t_z[] = { 1.0/6.0*PI, 5.0/6.0*PI, 9.0/6.0*PI, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case QUAD: {
		static const double t_e[] = { 0.0, 0.0, 0,0, 0.0, };
		static const double t_z[] = { PI,  0.0, 1.5*PI, 0.5*PI, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case TET: {
		static const double t_e[] = { THETA_E_TET-0.5*PI, THETA_E_TET-0.5*PI, THETA_E_TET-0.5*PI, 0.5*PI, };
		static const double t_z[] = { 1.0/6.0*PI, 5.0/6.0*PI, 9.0/6.0*PI, 0.0, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case HEX: {
		static const double t_e[] = { 0.0, 0.0, 0,0, 0.0, 0.5*PI, 1.5*PI, };
		static const double t_z[] = { PI,  0.0, 1.5*PI, 0.5*PI, 0.0, PI, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case WEDGE: {
		static const double t_e[] = { 0.0, 0.0, 0.0, 0.5*PI, 1.5*PI, };
		static const double t_z[] = { 1.0/6.0*PI, 5.0/6.0*PI, 9.0/6.0*PI, 0.0, 0.0, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} case PYR: {
		static const double t_e[] =
			{ THETA_E_PYR-0.5*PI, THETA_E_PYR-0.5*PI, THETA_E_PYR-0.5*PI, THETA_E_PYR-0.5*PI, 0.5*PI, };
		static const double t_z[] = { PI, 0.0, 1.5*PI, 0.5*PI, 0.0, };
		theta_eta  = t_e;
		theta_zeta = t_z;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.\n",e_type);
		break;
	}

	const int d   = e_info->d,
	          n_f = e_info->n_f;

	struct Multiarray_d* normals = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){n_f,d}); // returned

	double* data = normals->data;
	for (int ind = 0, f = 0; f < n_f; ++f) {
		data[ind++] = cos(theta_eta[f])*cos(theta_zeta[f]);
		if (d > 1)
			data[ind++] = cos(theta_eta[f])*sin(theta_zeta[f]);
		if (d > 2)
			data[ind++] = -sin(theta_eta[f]);
	}
	return normals;
}
