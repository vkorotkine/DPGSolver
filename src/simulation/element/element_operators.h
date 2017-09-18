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

struct Simulation;
struct const_Element;

/// Container for a Multiarray of \ref Cubature\* data.
struct Multiarray_Cubature {
	int order;          ///< Defined in \ref Multiarray_d.
	ptrdiff_t* extents; ///< Defined in \ref Multiarray_d.

	bool owns_data;              ///< Defined in \ref Multiarray_d.
	struct const_Cubature* data; ///< Defined in \ref Multiarray_d.
};

/// `const` version of \ref Multiarray_Cubature.
struct const_Multiarray_Cubature {
	const int order;               ///< Defined in \ref Multiarray_d.
	const ptrdiff_t*const extents; ///< Defined in \ref Multiarray_d.

	const bool owns_data;                   ///< Defined in \ref Multiarray_d.
	const struct const_Cubature*const data; ///< Defined in \ref Multiarray_d.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref const_Vector_i\* holding the extents for an operator to be subsequently constructed.
 *  \return Standard. */
const struct const_Vector_i* constructor_operator_extents_const_Vector_i
	(const struct Simulation* sim,        ///< \ref Simulation.
	 const struct const_Element* element, ///< \ref const_Element.
	 const int op_type                    ///< The operator type.
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
	(const struct Simulation* sim,          ///< \ref Simulation.
	 const struct const_Element* element,   ///< \ref const_Element.
	 const struct const_Vector_i* ext_v1_V, ///< The \ref Vector_i\* of operator extents.
	 const int cub_entity,                  ///< The entity for which the cubature nodes are used.
	 const int n_skip                       /**< The number of values to skip in `ext_v1_V`. Always skip the dimension
	                                         *   for differentiation operators as the nodes do not change. */
	);

/// \brief Destructor for a \ref const_Multiarray_Cubature\* container.
void destructor_const_Multiarray_Cubature
	(const struct const_Multiarray_Cubature*const a ///< Standard.
	);

#endif // DPG__element_operators_h__INCLUDED
