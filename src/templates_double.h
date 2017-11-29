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

#ifndef DPG__templates_double_h__INCLUDED
#define DPG__templates_double_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the `double` data type.
 */

#define TYPE_REAL 1 ///< Parameter used as a flag to enable functionality required for real functions.

///\{ \name Types.
#define Type double ///< Type parameter.

#define Multiarray_T       Multiarray_d       ///< Multiarray_\*.
#define const_Multiarray_T const_Multiarray_d ///< const_Multiarray_\*.
#define Multiarray_R       Multiarray_d       ///< Real Multiarray_\*.
#define const_Multiarray_R const_Multiarray_d ///< Real const_Multiarray_\*.
///\}

///\{ \name Boundary related parameters.
///\}

///\{ \name General functions.
#define constructor_empty_Multiarray_T      constructor_empty_Multiarray_d      ///< Standard.
#define constructor_copy_const_Multiarray_T constructor_copy_const_Multiarray_d ///< Standard.
#define destructor_const_Multiarray_T       destructor_const_Multiarray_d       ///< Standard.

#define set_to_value_Multiarray_T set_to_value_Multiarray_d ///< Standard.

#define get_col_const_Multiarray_T get_col_const_Multiarray_d ///< Standard.
#define get_col_Multiarray_T       get_col_Multiarray_d       ///< Standard.

#define abs_T fabs  ///< Standard.
#define real_T      ///< Standard.
#define sqrt_T sqrt ///< Standard.
#define pow_T  pow  ///< Standard.
///\}

#endif // DPG__templates_double_h__INCLUDED
