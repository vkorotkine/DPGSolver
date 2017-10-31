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
 *  Operator names take the general form: [type][0](7)_[1][2][3]_[4][5][6] where entries in [square brackets] are
 *  required and those in (round brackets) are optional:
 *  - type:
 *  	- cv: coefficients to values
 *  	- cc: coefficients to coefficients
 *  	- vv: values       to values
 *  	- vc: values       to coefficients
 *  	- tw: 't'est functions with 'w'eights.
 *  \note For collocated schemes, `tw` operators also include premultiplication by the inverse weights.
 *  - [0]:   value for the order of differentiation (0 for no differentiation).
 *  - [1/4]: character denoting the type of computational element (ce):
 *  	- 'v'olume
 *  	- 'f'ace
 *  	- 'e'dge
 *  - [2/5]: character denoting the kind of basis/nodes to be used (kind):
 *  	- 'g'eometry
 *  	- 'm'etric
 *  	- 'v'ertex (Only available as p1 nodes [5])
 *  	- 'p'lotting (Not available as a basis [2])
 *  	- 'c'ubature (Not available as a basis [2])
 *  	- 's'olution
 *  	- 'r': g'r'adient
 *  - (3/6): character denoting whether the basis/nodes is meant to be used for straight or curved elements (sc):
 *  	- 's'traight
 *  	- 'c'urved
 *  	- 'A'll
 *  - (7): if present, character specifying which basis was used for the operator:
 *  	- r: reference
 *
 *  The optional straight/curved parameters **must** be replaced with 'A'll if not present when passed to the
 *  constructor function.
 *
 *  Each operator also has an associated range with a maximum order of \ref OP_ORDER_MAX with the following parameters
 *  (d)(ce_o)(ce_i)[h_o][h_i][p_o][p_i], where entries in square and round brackets are once again required and
 *  optional, respectively.
 *
 *  When an operator is declared as an array of size 2, the two indices are used to denote the straight (0) and curved
 *  (1) operators, respectively.
 */

#include <stddef.h>
#include <stdbool.h>

struct Simulation;
struct const_Element;
struct Multiarray_Operator;
struct Operator_Info;

/** \brief Function pointer to a function which will set the operators.
 *
 *  \param ind_values The index of the first row of \ref Operator_Info::values_op to use.
 *  \param op         The multiarray of operators.
 *  \param op_info    \ref Operator_Info.
 *  \param sim        \ref Simulation.
 */
typedef void (*set_operator_fptr)
	(ptrdiff_t*const ind_values,
	 const struct Multiarray_Operator* op,
	 const struct Operator_Info* op_info,
	 const struct Simulation* sim
	);

/// Container specifying the effect ("from which input to which output") of the operator application.
struct Op_IO {
	const char ce,   ///< The computational element.
	           kind, ///< The kind of basis/nodes.
	           sc;   ///< Indication of straight/curved.

	const int ce_op, ///< The computational element index of the operator.
	          h_op,  ///< The h-refinement index of the operator.
	          p_op;  ///< The polynomial order index of the operator (**Not the order of the basis/cubature rule**).

	const int s_type; ///< \ref Element::s_type.
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

	/// The extents of the associated \ref Multiarray_Matrix_d\* of operators.
	const struct const_Vector_i* extents_op;

	const struct const_Matrix_i* values_op; ///< The values of d, f, h, p_in, and p_out for each operator.

	set_operator_fptr set_operator; ///< \ref set_operator_fptr.
};

// Interface functions ********************************************************************************************** //

/** \brief Constructor for a \ref Multiarray_Operator\* of operators.
 *  \return Standard. */
const struct Multiarray_Operator* constructor_operators
	(const char*const name_type,          ///< The name of the operator type (including differentiation index).
	 const char*const name_in,            ///< The name of the operator input.
	 const char*const name_out,           ///< The name of the operator output.
	 const char*const name_range,         ///< The name of the operator range.
	 const int p_ref[2],                  ///< Defined in \ref Operator_Info.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim         ///< \ref Simulation.
	);

/** \brief Constructor for a \ref const_Multiarray_Vector_i\* of 'n'ode 'c'orrespondence vectors.
 *  \return See brief. */
const struct const_Multiarray_Vector_i* constructor_operators_nc
	(const int ind_f_elem,                ///< Index of \ref Element::face_element.
	 const char*const name_in,            ///< Defined for \ref constructor_operators.
	 const char*const name_out,           ///< Defined for \ref constructor_operators.
	 const char*const name_range,         ///< Defined for \ref constructor_operators.
	 const int p_ref[2],                  ///< Defined for \ref constructor_operators.
	 const struct const_Element* element, ///< Defined for \ref constructor_operators.
	 const struct Simulation* sim         ///< Defined for \ref constructor_operators.
	);

/** \brief Constructor for a \ref const_Multiarray_Vector_d\* of cubature 'w'eight vectors.
 *  \return See brief. */
const struct const_Multiarray_Vector_d* constructor_operators_w
	(const char*const name_in,            ///< Defined for \ref constructor_operators.
	 const char*const name_out,           ///< Defined for \ref constructor_operators.
	 const char*const name_range,         ///< Defined for \ref constructor_operators.
	 const int p_ref[2],                  ///< Defined for \ref constructor_operators.
	 const struct const_Element* element, ///< Defined for \ref constructor_operators.
	 const struct Simulation* sim         ///< Defined for \ref constructor_operators.
	);

/** \brief Constructor for the \ref Operator_Info\* having the given inputs.
 *  \return Standard. */
struct Operator_Info* constructor_Operator_Info
	(const char*const name_type,         ///< Defined for \ref constructor_operators.
	 const char*const name_in,           ///< Defined for \ref constructor_operators.
	 const char*const name_out,          ///< Defined for \ref constructor_operators.
	 const char*const name_range,        ///< Defined for \ref constructor_operators.
	 const int p_ref[2],                 ///< Defined for \ref constructor_operators.
	 const struct const_Element* element ///< Defined for \ref constructor_operators.
	);

/// \brief Destructor for a \ref Operator_Info\*.
void destructor_Operator_Info
	(struct Operator_Info* op_ranges ///< Standard.
	);

/** \brief Compute the order of the basis based on the reference order and the kind of operator.
 *  \return See brief. */
int compute_p_basis
	(const struct Op_IO* op_io,   ///< \ref Op_IO.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Check if information can be lost while performing the operation.
 *  \return `true` if: `p_i > p_o` or `h_i > h_o`; `false` otherwise. */
bool check_op_info_loss
	(const int*const op_values ///< Values for the operator indices.
	);

/** \brief Constructor for a \ref Vector_i\* of indices for the current operator.
 *  \return See brief. */
const struct const_Vector_i* constructor_indices_Vector_i
	(const int ext_0_expected,     ///< The expected value of ext_0. May be set to '-1' if checking is not desired.
	 const int* op_values,         ///< The operator values.
	 const bool*const indices_skip ///< Indices to skip (if not NULL).
	);

#endif // DPG__element_operators_h__INCLUDED
