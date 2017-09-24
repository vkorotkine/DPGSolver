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
 *
 *  \todo Remove unused functions/classes.
 *
 *  Operator names take the general form: [type][0]_[1][2](3)_[4][5](6) where entries in [square brackets] are required
 *  and those in (round brackets) are optional:
 *  - type:
 *  	- cv: coefficients to values
 *  	- cc: coefficients to coefficients
 *  	- vv: values       to values
 *  	- vc: values       to coefficients
 *  - [0]:   value for the order of differentiation (0 for no differentiation).
 *  - [1/4]: character denoting the type of computational element (ce):
 *  	- v: volume
 *  	- f: face
 *  	- e: edge
 *  - [2/5]: character denoting the kind of basis/cubature to be used (kind):
 *  	- g: geometry
 *  	- m: metric
 *  	- s: solution
 *  	- v: vertex (Only available as p1 nodes [5])
 *  	- p: plotting (Not available as a basis [2])
 *  	- c: cubature (Not available as a basis [2])
 *  - (3/6): character denoting whether the basis/cubature is meant to be used for straight or curved elements (sc):
 *  	- s: straight
 *  	- c: curved
 *
 *  The optional straight/curved parameters **must** be replaced with 'A'll if not present when passed to the
 *  constructor function.
 *
 *  Each operator also has an associated range with a maximum order of \ref OP_ORDER_MAX with the following parameters
 *  (d)(f)[h_o][h_i][p_o][p_i], where entries in square and round brackets are once again required and optional,
 *  respectively.
 */

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct const_Element;

/// Container specifying the effect ("from which input to which output" ) of the operator application.
struct Op_IO {
	const char ce,   ///< The computational element.
	           kind, ///< The kind of basis/cubature.
	           sc;   ///< Indication of straight/curved.
};

/** Container for operator range related information.
 *  For available options for the parameters see the comments in \ref element_operators.h.
 */
struct Operator_Info {
	const struct const_Element* element; ///< \ref Element.

	const int op_type; ///< The type of operator.

	const struct Op_IO op_io[2]; ///< \ref Op_IO for each of input/output.

	const int range_d,  ///< Range of dimensions (For differentiation operators).
	          range_ce, ///< Range of computational elements.
	          range_h,  ///< Range of h-refinement related operators.
	          range_p;  ///< Range of polynomial orders.

	const int p_ref[2]; ///< Reference polynomial orders from \ref Simulation.

	const struct const_Vector_i* extents_cub; ///< The extents of the associated \ref Multiarray_Cubature\*.

	/// The extents of the associated \ref Multiarray_Matrix_d\* of operators.
	const struct const_Vector_i* extents_op;

	const struct const_Matrix_i* values_op ///< The values of d, f, h, p_in, and p_out for each operator.
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

/** \brief Constructor for a \ref const_Multiarray_Matrix_d\* of operators.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_operators
	(const char*const name_type,          ///< The name of the operator type (including differentiation index).
	 const char*const name_in,            ///< The name of the operator input.
	 const char*const name_out,           ///< The name of the operator output.
	 const char*const name_range,         ///< The name of the operator range.
	 const int p_ref[2],                  ///< Defined in \ref Operator_Info.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

#endif // DPG__element_operators_h__INCLUDED
