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

#ifndef DPG__element_solver_dpg_h__INCLUDED
#define DPG__element_solver_dpg_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref DPG_Solver_Element container and associated functions.
 */

#include "element_solver.h"

struct Simulation;

/// \brief Container for data relating to the dpg solver element.
struct DPG_Solver_Element {
	const struct Solver_Element element; ///< Base \ref Solver_Element.

	// Volume rlhs
	const struct Multiarray_Operator* cv0_vt_vc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vt_vc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cvcv0_vs_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cvcv1_vt_vc[2]; ///< See notation in \ref element_operators.h.

	// Face rlhs
	const struct Multiarray_Operator* cv0_ff_fc[2]; ///< See notation in \ref element_operators.h.

	// Enforcing conservation
	const struct const_Multiarray_Vector_d* ones_coef_vt; ///< \See \ref constructor_operators_ones_coef.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref DPG_Solver_Element.
void constructor_derived_DPG_Solver_Element
	(struct Element* element,     ///< \ref DPG_Solver_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref DPG_Solver_Element.
void destructor_derived_DPG_Solver_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solver_dpg_h__INCLUDED
