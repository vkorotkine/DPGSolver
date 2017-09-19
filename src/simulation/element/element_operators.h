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

#ifndef DPG__element_operators_h__INCLUDED
#define DPG__element_operators_h__INCLUDED
/** \file
 *  \brief Provides the functions relating to element operators.
 */

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct const_Element;

/// Container for operator range related information.
struct Operator_Info {
	/** The type of operator. Options:
	 *  - 'T'ransform:   coefficients to coefficients
	 *  - 'E'valuate:    coefficients to values
	 *  - 'I'nterpolate: values       to values
	 *  - 'P'roject:     values       to coefficients
	 */
	char op_type;
// Potentially remove.

	const int range_d, ///< Range of dimensions (For differentiation operators).
	          range_f, ///< Range of faces.
	          range_p, ///< Range of orders.
	          range_h; ///< Range of h-refinement related operators.

	const int cub_type; ///< The type of cubature.

	const int p_ref[2]; ///< Reference polynomial orders from \ref Simulation.
};

/// Container for a Multiarray of \ref Cubature\* data.
struct Multiarray_Cubature {
	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data;                    ///< Defined in \ref Multiarray_d.
	const struct const_Cubature** data; ///< Defined in \ref Multiarray_d.
};

/// `const` version of \ref Multiarray_Cubature.
struct const_Multiarray_Cubature {
	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;                         ///< Defined in \ref Multiarray_d.
	const struct const_Cubature*const*const data; ///< Defined in \ref Multiarray_d.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for the \ref Operator_Info\* having the given inputs.
 *  \return Standard. */
struct Operator_Info* constructor_Operator_Info
	(const int range_d,  ///< Defined in \ref Operator_Info.
	 const int range_f,  ///< Defined in \ref Operator_Info.
	 const int range_p,  ///< Defined in \ref Operator_Info.
	 const int range_h,  ///< Defined in \ref Operator_Info.
	 const int cub_type, ///< Defined in \ref Operator_Info.
	 const int p_ref[2]  ///< Defined in \ref Operator_Info.
	);

/// \brief Destructor for a \ref Operator_Info\*.
void destructor_Operator_Info
	(struct Operator_Info* op_ranges ///< Standard.
	);

/** \brief Constructor for a \ref const_Multiarray_Cubature\* holding the cubature nodes (and weights if applicable) for
 *         the range of supported hp adaptive operators.
 *  \return Standard.
 *
 *  The required nodes are determined from the last `n_hp` values in the `ext_v1_V` vector. Cubature nodes on
 *  element sub-regions are computed by:
 *  1. Computing the barycentric coordinates of the nodes on the complete reference element;
 *  2. Multiplying with the appropriate vertices of the element sub-regions.
 */
const struct const_Multiarray_Cubature* constructor_const_Multiarray_Cubature
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Operator_Info* op_info  ///< \ref Operator_Info.
	);

/// \brief Destructor for a \ref const_Multiarray_Cubature\* container.
void destructor_const_Multiarray_Cubature
	(const struct const_Multiarray_Cubature*const a ///< Standard.
	);

#endif // DPG__element_operators_h__INCLUDED
