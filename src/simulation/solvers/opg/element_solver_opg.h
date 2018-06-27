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

#ifndef DPG__element_solver_opg_h__INCLUDED
#define DPG__element_solver_opg_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref OPG_Solver_Element container and associated functions.
 */

#include "element_solver.h"

struct Simulation;

/// \brief Container for data relating to the opg solver element.
struct OPG_Solver_Element {
	const struct Solver_Element element; ///< Base \ref Solver_Element.

	// Volume
	const struct Multiarray_Operator* cv0_vt_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vt_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vc0_vs_vs;    ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vt_vs;    ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vg_vt[2]; ///< See notation in \ref element_operators.h.

	// Face
	const struct Multiarray_Operator* cv0_vt_fc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vt_fc[2]; ///< See notation in \ref element_operators.h.

	// Tensor-product sub-operators.
	const struct Multiarray_Operator* cv0_vt_vs;    ///< See notation in \ref element_operators.h.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref OPG_Solver_Element.
void constructor_derived_OPG_Solver_Element
	(struct Element* element,     ///< \ref OPG_Solver_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref OPG_Solver_Element.
void destructor_derived_OPG_Solver_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_solver_opg_h__INCLUDED
