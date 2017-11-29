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

#ifndef DPG__templates_double_complex_h__INCLUDED
#define DPG__templates_double_complex_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the `double complex` data type.
 *
 *  See \ref templates_double.h for description of parameters.
 */

#define TYPE_COMPLEX 1 ///< See brief.

///\{ \name Types.
#define Type double complex ///< See brief.

#define Multiarray_T       Multiarray_c       ///< See brief.
#define const_Multiarray_T const_Multiarray_c ///< See brief.
#define Multiarray_R       Multiarray_d       ///< See brief.
#define const_Multiarray_R const_Multiarray_d ///< See brief.
///\}

///\{ \name General functions.
#define constructor_empty_Multiarray_T      constructor_empty_Multiarray_c      ///< See brief.
#define constructor_copy_const_Multiarray_T constructor_copy_const_Multiarray_c ///< See brief.
#define destructor_const_Multiarray_T       destructor_const_Multiarray_c       ///< See brief.

#define set_to_value_Multiarray_T           set_to_value_Multiarray_c ///< See brief.

#define get_col_const_Multiarray_T get_col_const_Multiarray_c ///< See brief.
#define get_col_Multiarray_T       get_col_Multiarray_c       ///< See brief.

#define abs_T  cabs  ///< See brief.
#define real_T creal ///< See brief.
#define sqrt_T csqrt ///< See brief.
#define pow_T  cpow  ///< See brief.
///\}

#endif // DPG__templates_double_complex_h__INCLUDED
