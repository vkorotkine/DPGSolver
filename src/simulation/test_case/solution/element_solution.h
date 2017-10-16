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

#ifndef DPG__element_solution_h__INCLUDED
#define DPG__element_solution_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Solution_Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref Solution_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the solution elements.
struct Solution_Element {
	struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* cv0_vgs_vss; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vsc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgs_vrs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vrc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vgs_fss; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_fsc; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgs_frs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_frc; ///< See notation in \ref element_operators.h.
};

/// \brief `const` version of the \ref Solution_Element container.
struct const_Solution_Element {
	const struct const_Element element; ///< Defined in \ref Solution_Element.

	const struct Multiarray_Operator*const cv0_vgs_vss; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgc_vsc; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgs_vrs; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgc_vrc; ///< Defined in \ref Solution_Element.

	const struct Multiarray_Operator*const cv0_vgs_fss; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgc_fsc; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgs_frs; ///< Defined in \ref Solution_Element.
	const struct Multiarray_Operator*const cv0_vgc_frc; ///< Defined in \ref Solution_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Solution_Element.
void constructor_derived_Solution_Element
	(struct Element* element,     ///< \ref Solution_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Solution_Element.
void destructor_derived_Solution_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solution_h__INCLUDED
