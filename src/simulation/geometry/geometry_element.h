// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__geometry_element_h__INCLUDED
#define DPG__geometry_element_h__INCLUDED
/**	\file
 *	\brief Provides the interface for the derived \ref Geometry_Element container and associated functions.
 *
 *	\note `const` and non-`const` versions of \ref Geometry_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the geometry elements.
struct Geometry_Element {
	struct const_Element element; ///< Base \ref const_Element.

};

/// \brief `const` version of the \ref Geometry_Element container.
struct const_Geometry_Element {
	struct const_Element element; ///< Base \ref const_Element.

};

// Constructor/Destructor functions ********************************************************************************* //

/** \brief Constructs the \ref Geometry_Element\*s.
 *	\return Standard. */
struct const_Intrusive_List* constructor_Geometry_Elements
	(struct Simulation*const sim ///< The \ref Simulation.
	);

/// \brief Destructs the \ref Geometry_Element\*s.
void destructor_Geometry_Elements
	(struct Intrusive_List* geometry_elements ///< Standard.
	);


#endif // DPG__geometry_element_h__INCLUDED
