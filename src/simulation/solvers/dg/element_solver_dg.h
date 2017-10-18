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

#ifndef DPG__element_solver_dg_h__INCLUDED
#define DPG__element_solver_dg_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref DG_Solver_Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref DG_Solver_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the dg solver element.
struct DG_Solver_Element {
	struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* cv0_vs_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vs_vcc; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* tw1_vs_vcs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* tw1_vs_vcc; ///< See notation in \ref element_operators.h.
};

/// \brief `const` version of the \ref DG_Solver_Element container.
struct const_DG_Solver_Element {
	const struct const_Element element; ///< Defined in \ref DG_Solver_Element.

	const struct Multiarray_Operator*const cv0_vs_vcs; ///< Defined in \ref DG_Solver_Element.
	const struct Multiarray_Operator*const cv0_vs_vcc; ///< Defined in \ref DG_Solver_Element.

	const struct Multiarray_Operator*const tw1_vs_vcs; ///< Defined in \ref DG_Solver_Element.
	const struct Multiarray_Operator*const tw1_vs_vcc; ///< Defined in \ref DG_Solver_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref DG_Solver_Element.
void constructor_derived_DG_Solver_Element
	(struct Element* element,     ///< \ref DG_Solver_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref DG_Solver_Element.
void destructor_derived_DG_Solver_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solver_dg_h__INCLUDED
