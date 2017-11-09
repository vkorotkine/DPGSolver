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

#ifndef DPG__element_solver_h__INCLUDED
#define DPG__element_solver_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Solver_Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref Solver_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the dg solver element.
struct Solver_Element {
	const struct const_Element element; ///< Base \ref const_Element.

	// Volume rlhs
	const struct Multiarray_Operator* cv0_vs_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* tw1_vt_vc[2]; ///< See notation in \ref element_operators.h.

	// Tensor-Product
	const struct Multiarray_Operator* tw0_vt_vc[2]; ///< See notation in \ref element_operators.h.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Solver_Element.
void constructor_derived_Solver_Element
	(struct Element* element,     ///< \ref Solver_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Solver_Element.
void destructor_derived_Solver_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solver_h__INCLUDED
