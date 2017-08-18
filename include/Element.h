// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Element_h__INCLUDED
#define DPG__Element_h__INCLUDED
/**	\file
 *	\todo Provides the interface for the base \ref Element container and associated functions.
 */

#include "Intrusive.h"

/// \brief Container for data relating to the base Elements.
struct Element {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const unsigned int type, ///< The element type.
	                   d,    ///< The dimension.
	                   n_f;  ///< The number of faces.

	const struct const_Multiarray_Vector_ui*const f_ve; ///< The correspondence between the (f)aces and (ve)rtices.

};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Element \ref Intrusive_List.
struct Intrusive_List* constructor_Element_List
	(const unsigned int d ///< The dimension.
	);

/// \brief Destructs the base \ref Element \ref Intrusive_List.
void destructor_Elements
	(struct Intrusive_List* Elements ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief See return.
 *	\return Pointer to a base \ref Element of the input `type`.
 */
struct Element* get_element_by_type
	(const struct Intrusive_List* elements, ///< The list of elements.
	 const unsigned int type                ///< Defined in \ref Element.
	);


#endif // DPG__Element_h__INCLUDED
