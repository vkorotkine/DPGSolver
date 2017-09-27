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

#ifndef DPG__element_operators_tp_h__INCLUDED
#define DPG__element_operators_tp_h__INCLUDED
/** \file
 *  \brief Provides the functions relating to tensor-product element operators.
 *
 *  See \ref element_operators.h for relevant comments.
 *
 *  \warning The limited number of input operators currently does not allow for certain tensor-product operators
 *           to be constructed.
 *  For example, a cv1_vgc_fcc WEDGE operator for h-adaptation could not be constructed as it would require 8 inputs:
 *  cv0_vgc_vcc, cv1_vgc_vcc, cv0_vgc_fcc, cv1_vgc_fcc operators for both LINE and TRI elements.
 *
 *  While this is not currently an issue as no such operator is used, it would perhaps be convenient to store the
 *  non-differentiation and differentiation operators as part of the same Multiarray, leading once again to 4 operators.
 */

struct Simulation;
struct const_Element;

/// Container for the sub-element operators used to construct the tensor-product operators.
struct Operators_TP {
	/** Sub-element operators.
	 *
	 *  The indices are over the sub-elements and sub-element operators respectively. For example:
	 *  	wedge_op[0][1] : tri  operator 1.
	 *  	wedge_op[1][0] : line operator 0.
	 */
	const struct const_Multiarray_Matrix_d* op[2][2];
};

// Interface functions ********************************************************************************************** //

/// \brief Set the entries of the input \ref Operators_TP container.
void set_operators_tp
	(struct Operators_TP* ops_tp,             ///< \ref Operators_TP.
	 const struct const_Multiarray_Matrix_d* op_00, ///< Sets \ref Operators_TP::op, index [0][0].
	 const struct const_Multiarray_Matrix_d* op_01, ///< Sets \ref Operators_TP::op, index [0][1].
	 const struct const_Multiarray_Matrix_d* op_10, ///< Sets \ref Operators_TP::op, index [1][0].
	 const struct const_Multiarray_Matrix_d* op_11  ///< Sets \ref Operators_TP::op, index [1][1].
	);

/** \brief Constructor for a \ref const_Multiarray_Matrix_d\* of tensor-product operators from sub-element operators.
 *  \return Standard. */
const struct const_Multiarray_Matrix_d* constructor_operators_tp
	(const char*const name_type,          ///< The name of the operator type (including differentiation index).
	 const char*const name_in,            ///< The name of the operator input.
	 const char*const name_out,           ///< The name of the operator output.
	 const char*const name_range,         ///< The name of the operator range.
	 const struct const_Element* element, ///< \ref const_Element.
	 const struct Simulation* sim,        ///< \ref Simulation.
	 const struct Operators_TP* ops_tp    ///< \ref Operators_TP.
	);

#endif // DPG__element_operators_tp_h__INCLUDED
