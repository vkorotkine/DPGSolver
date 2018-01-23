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
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the solution elements.
struct Solution_Element {
	const struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* cv0_vg_vs[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_vr[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vs_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vs_vs;    ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vc0_vs_vs;    ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vc0_vr_vr;    ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vg_ff[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vm_ff[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_fc[2]; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* vc0_ff_ff; ///< See notation in \ref element_operators.h.

	// Error computation
	const struct const_Multiarray_Vector_d* w_vc[2]; ///< Weights for 'v'olume 'c'ubature 's'traight/'c'urved.

	// Tensor-product sub-operators.
	const struct Multiarray_Operator* cv0_vg_vf[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* vv0_vm_vf[2]; ///< See notation in \ref element_operators.h.
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
