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

	const struct Multiarray_Operator* cv1_vgs_vcs; ///< See notation in \todo [Ref here].
	const struct Multiarray_Operator* cv1_vgc_vcc; ///< See notation in \todo [Ref here].
};

/// \brief `const` version of the \ref Geometry_Element container.
struct const_Geometry_Element {
	const struct const_Element element; ///< Defined in \ref Geometry_Element.

	const struct Multiarray_Operator*const cv1_vgs_vcs; ///< Defined in \ref Geometry_Element.
	const struct Multiarray_Operator*const cv1_vgc_vcc; ///< Defined in \ref Geometry_Element.
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
