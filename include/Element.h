// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Element_h__INCLUDED
#define DPG__Element_h__INCLUDED
/**	\file
 *	Provides the interface for the base \ref Element container and associated functions.
 *
 *	\note `const` and non-`const` versions of \ref Element must have identical members and layout.
 */

#include "Intrusive.h"

/// \brief Container for data relating to the base Elements.
struct Element {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const int type, ///< The element type.
	          d,    ///< The dimension.
	          n_f;  ///< The number of faces.

	const struct const_Multiarray_Vector_i*const f_ve; ///< The correspondence between the (f)aces and (ve)rtices.

};

/// \brief `const` version of the base Element container.
struct const_Element {
	const struct const_Intrusive_Link lnk;

	const int type,
	          d,
	          n_f;

	const struct const_Multiarray_Vector_i*const f_ve;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Element \ref Intrusive_List.
struct const_Intrusive_List* constructor_Element_List
	(const int d ///< The dimension.
	);

/// \brief Destructs the base \ref Element \ref Intrusive_List.
void destructor_Elements
	(struct Intrusive_List* Elements ///< Standard.
	);

/// \brief Cast from \ref const_Element\* to `const` \ref const_Element `*const`.
void const_cast_const_Element
	(const struct const_Element*const* dest, ///< Destination.
	 const struct const_Element*const src    ///< Source.
	);

// Helper functions ************************************************************************************************* //

/** \brief See return.
 *	\return Pointer to a base \ref Element of the input `type`.
 */
struct const_Element* get_element_by_type
	(const struct const_Intrusive_List*const elements, ///< The list of elements.
	 const int type                                    ///< Defined in \ref Element.
	);

/** \brief See return.
 *	\return Pointer to a base \ref Element corresponding to the specified face of the input volume.
 */
struct const_Element* get_element_by_face
	(const struct const_Element*const element, ///< The element corresponding to the volume.
	 const int lf                              ///< The index of the local face.
	);

#endif // DPG__Element_h__INCLUDED
