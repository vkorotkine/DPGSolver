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

#ifndef DPG__element_solver_fr_split_form_h__INCLUDED
#define DPG__element_solver_fr_split_form_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref FR_SPLIT_FORM_Solver_Element container and associated functions.
 */

#include "element_solver.h"
struct Simulation;

/// \brief Container for data relating to the fr_split_form solver element.
struct FRSF_Solver_Element {
	const struct Solver_Element element; ///< Base \ref Solver_Element.

	const struct const_Matrix_d* mass;
	const struct const_Matrix_d* D_frsf;

	const struct const_Matrix_d* cv0_vs_vc;
	const struct const_Matrix_d* cv1_vs_vc;
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref FRSF_Solver_Element.
void constructor_derived_FRSF_Solver_Element
	(struct Element* element,     ///< \ref FRSF_Solver_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref FRSF_Solver_Element.
void destructor_derived_FRSF_Solver_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solver_fr_split_form_h__INCLUDED
