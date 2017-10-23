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

#ifndef DPG__element_h__INCLUDED
#define DPG__element_h__INCLUDED
/** \file
 *  \brief Provides the interface for the base \ref Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref Element must have identical members and layout.
 */

#include <stdbool.h>
#include "intrusive.h"

struct const_Vector_i;

/// \brief Container for data relating to the base Elements.
struct Element {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	bool present; ///< Whether the element type is or may become present in the simulation.

	int type,   ///< The element type.
	    s_type, ///< The element super type. Options: Tensor-Product, SImplex, PYRamid, WEDGE.
	    d,      ///< The dimension.
	    n_ve,   ///< The number of vertices.
	    n_e,    ///< The number of edges.
	    n_f;    ///< The number of faces.

	int n_ref_max_v, ///< Maximum number of volume h-refinements.
	    n_ref_max_f, ///< Maximum number of face h-refinements.
	    n_ref_max_e; ///< Maximum number of edge h-refinements.

	struct Multiarray_Vector_i* f_ve; ///< The correspondence between the (f)aces and (ve)rtices.

	struct Multiarray_d* normals; ///< The reference unit normal vectors for each of the element faces.

	struct Element* sub_element[2];  ///< Sub-element(s) used to form the tensor-product element (if applicable).
	struct Element* face_element[2]; ///< Element(s) associated with the faces of the \ref Element.

	void* derived; ///< Pointer to the element currently derived from the base element.
};

/// \brief `const` version of the \ref Element container.
struct const_Element {
	const struct const_Intrusive_Link lnk; ///< Defined in \ref Element.

	const bool present; ///< Defined in \ref Element.

	const int type,   ///< Defined in \ref Element.
	          s_type, ///< Defined in \ref Element.
	          d,      ///< Defined in \ref Element.
	          n_ve,   ///< Defined in \ref Element.
	          n_e,    ///< Defined in \ref Element.
	          n_f;    ///< Defined in \ref Element.

	const int n_ref_max_v, ///< Defined in \ref Element.
	          n_ref_max_f, ///< Defined in \ref Element.
	          n_ref_max_e; ///< Defined in \ref Element.

	const struct const_Multiarray_Vector_i*const f_ve; ///< Defined in \ref Element.

	const struct const_Multiarray_d* normals; ///< Defined in \ref Element.

	const struct const_Element*const sub_element[2];  ///< Defined in \ref Element.
	const struct const_Element*const face_element[2]; ///< Defined in \ref Element.

	const void*const derived; ///< Defined in \ref Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs the \ref Element\*s.
 *  \return Standard. */
struct const_Intrusive_List* constructor_Elements
	(const int d ///< The dimension.
	);

/// \brief Destructs the \ref Element\*s.
void destructor_Elements
	(struct Intrusive_List* elements ///< Standard.
	);

/// \brief `const` version of \ref destructor_Elements.
void destructor_const_Elements
	(const struct const_Intrusive_List* elements ///< Standard.
	);

/// \brief Destructor for an individual \ref Element.
void destructor_Element
	(struct Element* element ///< Standard.
	);

/// \brief Cast from \ref const_Element\* to `const` \ref const_Element `*const`.
void const_cast_const_Element
	(const struct const_Element*const* dest, ///< Destination.
	 const struct const_Element*const src    ///< Source.
	);

/// \brief Set the pointers to \ref Element::sub_element when applicable (QUAD, HEX, WEDGE).
void set_tp_sub_elements
	(struct Intrusive_List* elements ///< The list of elements.
	);

/// \brief Set the pointers to \ref Element::face_element.
void set_face_elements
	(struct Intrusive_List* elements ///< The list of elements.
	);

// Helper functions ************************************************************************************************* //

/// \brief Set \ref Element::present to `true` for element types which are present in the \ref Mesh_Data::elem_types.
void set_elements_present
	(const struct const_Intrusive_List* elements, ///< \ref Simulation::elements.
	 const struct const_Vector_i*const elem_types ///< \ref Mesh_Data::elem_types.
	);

/// \brief Remove \ref Element\*s which are not present from the list, except for the POINT element.
void remove_absent_Elements
	(const struct const_Intrusive_List* elements ///< \ref Simulation::elements.
	);

/** \brief See return.
 *  \return Pointer to a \ref Element of the input `type`. */
const struct const_Element* get_element_by_type
	(const struct const_Intrusive_List*const elements, ///< The list of elements.
	 const int type                                    ///< Defined in \ref Element.
	);

/** \brief See return.
 *  \return Pointer to a \ref Element corresponding to the specified face of the input volume. */
const struct const_Element* get_element_by_face
	(const struct const_Element*const element, ///< The element corresponding to the volume.
	 const int lf                              ///< The index of the local face.
	);

/** \brief Check whether wedges are present in the list of elements.
 *  \return `true` if yes, `false` otherwise. */
bool wedges_present
	(const struct const_Intrusive_List*const elements ///< The list of elements.
	);

/** \brief Compute the element type for the given input parameters.
 *  \return See brief. */
int compute_elem_type_sub_ce
	(const int e_type, ///< \ref Element::type.
	 const char ce,    ///< The type of computational element.
	 const int ind_ce  ///< The index of the computational element.
	);

/** \brief Compute the super type for the input element type.
 *  \return See brief. */
int compute_super_from_elem_type
	(const int e_type ///< \ref Element::type.
	);

/** \brief Compute the element type from the input super type and dimension.
 *  \return See brief. */
int compute_elem_from_super_type
	(const int s_type, ///< \ref Element::s_type.
	 const int d       ///< \ref Element::d.
	);

#endif // DPG__element_h__INCLUDED
