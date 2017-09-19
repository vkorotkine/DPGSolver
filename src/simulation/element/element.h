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

#include "intrusive.h"

/// \brief Container for data relating to the base Elements.
struct Element {
/// \todo Remove `const`s here.
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const int type,   ///< The element type.
	          s_type, ///< The element super type. Options: Tensor-Product, SImplex, PYRamid, WEDGE.
	          d,      ///< The dimension.
	          n_ve,   ///< The number of vertices.
	          n_f;    ///< The number of faces.

	const int n_ref_max; ///< Maximum number of h-refinements.

	const struct const_Multiarray_Vector_i*const f_ve; ///< The correspondence between the (f)aces and (ve)rtices.
};

/// \brief `const` version of the \ref Element container.
struct const_Element {
	const struct const_Intrusive_Link lnk; ///< Defined in \ref Element.

	const int type,   ///< Defined in \ref Element.
	          s_type, ///< Defined in \ref Element.
	          d,      ///< Defined in \ref Element.
	          n_ve,   ///< Defined in \ref Element.
	          n_f;    ///< Defined in \ref Element.

	const int n_ref_max; ///< Defined in \ref Element.

	const struct const_Multiarray_Vector_i*const f_ve; ///< Defined in \ref Element.
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

// Helper functions ************************************************************************************************* //

/** \brief See return.
 *  \return Pointer to a \ref Element of the input `type`. */
struct const_Element* get_element_by_type
	(const struct const_Intrusive_List*const elements, ///< The list of elements.
	 const int type                                    ///< Defined in \ref Element.
	);

/** \brief See return.
 *  \return Pointer to a \ref Element corresponding to the specified face of the input volume.
 */
struct const_Element* get_element_by_face
	(const struct const_Element*const element, ///< The element corresponding to the volume.
	 const int lf                              ///< The index of the local face.
	);

#endif // DPG__element_h__INCLUDED
