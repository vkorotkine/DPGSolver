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

#ifndef DPG__nodes_plotting_h__INCLUDED
#define DPG__nodes_plotting_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions computing the reference element coordinates for high-order plotting.
 */

#include <stdbool.h>
#include "nodes.h"

/// Container for reference element coordinates and cubature related information.
struct Plotting_Nodes {
	struct Nodes nodes; ///< \ref Nodes.

	/// The multiarray of node connectivity for the sub-elements within the element.
	struct Multiarray_Vector_i* connect;

	/// The multiarray of node connectivity for the sub-element edges of the element.
	struct Multiarray_Vector_i* connect_e;

	/// Indices of the linear VTK element types used by paraview. See the Figure 2 in $ROOT/doc/VTK_file_formats.pdf.
	struct Vector_i* vtk_types;
};

/// `const` version of \ref Plotting_Nodes.
struct const_Plotting_Nodes {
	const struct const_Nodes nodes; ///< Defined in \ref Plotting_Nodes.

	const struct const_Multiarray_Vector_i* connect;   ///< Defined in \ref Plotting_Nodes.
	const struct const_Multiarray_Vector_i* connect_e; ///< Defined in \ref Plotting_Nodes.

	const struct const_Vector_i* vtk_types; ///< Defined in \ref Plotting_Nodes.
};

// Constructor functions ******************************************************************************************** //

/** \brief Constructor for a \ref const_Plotting_Nodes container for the given order and element type.
 *  \return See brief. */
const struct const_Plotting_Nodes* constructor_const_Plotting_Nodes
	(const int p,     ///< Defined in \ref constructor_Nodes_fptr.
	 const int e_type ///< \ref Element::type.
	);

/// \brief Destructor for a \ref const_Plotting_Nodes\* container.
void destructor_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const p_nodes ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/// \brief Print a \ref Plotting_Nodes container to the terminal displaying entries below the tolerance as 0.0.
void print_Plotting_Nodes_tol
	(const struct Plotting_Nodes*const p_nodes, ///< Standard.
	 const double tol                           ///< The tolerance.
	);

/// \brief `const` version of \ref print_Plotting_Nodes_tol.
void print_const_Plotting_Nodes_tol
	(const struct const_Plotting_Nodes*const p_nodes, ///< Defined for \ref print_Plotting_Nodes_tol.
	 const double tol                                 ///< Defined for \ref print_Plotting_Nodes_tol.
	);

/** \brief Print a \ref Plotting_Nodes container to the terminal calling \ref print_Plotting_Nodes_tol with a default
 *         tolerance. */
void print_Plotting_Nodes
	(const struct Plotting_Nodes*const p_nodes ///< Standard.
	);

/// \brief `const` version of \ref print_Plotting_Nodes.
void print_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const p_nodes ///< Standard.
	);

#endif // DPG__nodes_plotting_h__INCLUDED
