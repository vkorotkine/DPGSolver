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

#ifndef DPG__element_plotting_h__INCLUDED
#define DPG__element_plotting_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Plotting_Element container and associated functions.
 *
 *  \note `const` and non-`const` versions of \ref Plotting_Element must have identical members and layout.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the plotting elements.
struct Plotting_Element {
	struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* cv0_vgs_vp; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vgc_vp; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cv0_vs_vp;  ///< See notation in \ref element_operators.h.

	int n_p;                                     ///< The number of stored plotting nodes.
	const struct const_Plotting_Nodes** p_nodes; ///< The array of \ref Plotting_Nodes for each order.
};

/// \brief `const` version of the \ref Plotting_Element container.
struct const_Plotting_Element {
	const struct const_Element element; ///< Defined in \ref Plotting_Element.

	const struct Multiarray_Operator*const cv0_vgs_vp; ///< Defined in \ref Plotting_Element.
	const struct Multiarray_Operator*const cv0_vgc_vp; ///< Defined in \ref Plotting_Element.
	const struct Multiarray_Operator*const cv0_vs_vp;  ///< Defined in \ref Plotting_Element.

	const int n_p;                                         ///< Defined in \ref Plotting_Element.
	const struct const_Plotting_Nodes*const*const p_nodes; ///< Defined in \ref Plotting_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Plotting_Element.
void constructor_derived_Plotting_Element
	(struct Element* element,     ///< \ref Plotting_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Plotting_Element.
void destructor_derived_Plotting_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_plotting_h__INCLUDED
