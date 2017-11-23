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

#ifndef DPG__element_error_h__INCLUDED
#define DPG__element_error_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Error_Element container and associated functions.
 */

#include "element_solution.h"

struct Simulation;

/// \brief Container for data relating to the error element.
struct Error_Element {
	const struct Solution_Element element; ///< Base \ref Solution_Element.

	const struct const_Multiarray_Vector_d* w_vc[2]; ///< Weights for 'v'olume 'c'ubature 's'traight/'c'urved.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Error_Element.
void constructor_derived_Error_Element
	(struct Element* element,     ///< \ref Error_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Error_Element.
void destructor_derived_Error_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_error_h__INCLUDED
