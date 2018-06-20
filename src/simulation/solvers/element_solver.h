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
 */

#include "element.h"
#include "element_adaptation.h"
#include "element_geometry.h"
#include "element_plotting.h"
#include "element_solution.h"

struct Simulation;

/// \brief Container for data relating to the solver element.
struct Solver_Element {
	const struct const_Element element; ///< Base \ref const_Element.

	const struct Adaptation_Element a_e; ///< \ref Adaptation_Element.
	const struct Geometry_Element   g_e; ///< \ref Geometry_Element.
	const struct Plotting_Element   p_e; ///< \ref Plotting_Element.
	const struct Solution_Element   s_e; ///< \ref Solution_Element.

	// Volume rlhs
	const struct Multiarray_Operator* cv0_vs_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vr_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* tw1_vt_vc[2]; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* cv0_vs_vs;    ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vr_vs;    ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_vs[2]; ///< See notation in \ref element_operators.h.

	const struct const_Multiarray_Vector_d* w_vc[2]; ///< Weights for 'v'olume 'c'ubature.

	const struct Multiarray_Operator* cv0_vt_vc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv1_vt_vc[2];   ///< See notation in \ref element_operators.h.

	// Face rlhs
	const struct Multiarray_Operator* cv0_vs_fc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vr_fc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* tw0_vt_fc[2];   ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_ff_fc[2];   ///< See notation in \ref element_operators.h.
	const struct const_Multiarray_Vector_i* nc_fc[2]; ///< Node correspondence face 'f'ace 'c'ubature.

	const struct const_Multiarray_Vector_d* w_fc[2]; ///< Weights for 'f'ace 'c'ubature.

	// Source rhs
	const struct Multiarray_Operator* cv0_vg_vc[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* tw0_vt_vc[2]; ///< See notation in \ref element_operators.h.

	// Positivity preservation
	const struct Multiarray_Operator* ccSB0_vs_vs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* ccBS0_vs_vs; ///< See notation in \ref element_operators.h.

	// CFL ramping
	const struct Multiarray_Operator* cv0_vg_vv[2]; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vg_ev[2]; ///< See notation in \ref element_operators.h.
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
