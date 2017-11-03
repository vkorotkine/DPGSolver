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

#ifndef DPG__test_support_operators_h__INCLUDED
#define DPG__test_support_operators_h__INCLUDED
/** \file
 *  \brief Provides `complex` version of functions related to the \ref Operator\* container.
 */

#include <stddef.h>

struct Operator;
struct Multiarray_c;

/** \brief `complex` version of \ref constructor_mm_NN1_Operator_const_Multiarray_d.
 *  \return See brief. */
const struct const_Multiarray_c* constructor_mm_NN1_Operator_const_Multiarray_c
	(const struct Operator* op,          ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	 const struct const_Multiarray_c* b, ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	 const char layout_c,                ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	 const char op_format,               ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	 const int order_sub_ma,             ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	 const ptrdiff_t* sub_inds_b         ///< Defined for \ref constructor_mm_NN1_Operator_const_Multiarray_d.
	);

/// \brief `complex` version of \ref mm_NNC_Operator_Multiarray_d.
void mm_NNC_Operator_Multiarray_c
	(const double alpha,                 ///< See brief.
	 const double beta,                  ///< See brief.
	 const struct Operator* op,          ///< See brief.
	 const struct const_Multiarray_c* b, ///< See brief.
	 struct Multiarray_c* c,             ///< See brief.
	 const char op_format,               ///< See brief.
	 const int order_sub_ma,             ///< See brief.
	 const ptrdiff_t* sub_inds_b,        ///< See brief.
	 const ptrdiff_t* sub_inds_c         ///< See brief.
	);

#endif // DPG__test_support_operators_h__INCLUDED
