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

#ifndef DPG__templates_multiarray_c_h__INCLUDED
#define DPG__templates_multiarray_c_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the `double complex` multiarray
 *         containers/functions.
 */

#include "templates_multiarray_constructors_c.h"

#define Multiarray_T       Multiarray_c       ///< Multiarray_\*.
#define Multiarray_R       Multiarray_d       ///< Real Multiarray_\*.
#define const_Multiarray_T const_Multiarray_c ///< const_Multiarray_\*.
#define const_Multiarray_R const_Multiarray_d ///< Real const_Multiarray_\*.


#define set_to_value_Multiarray_T set_to_value_Multiarray_c ///< Standard.

#define get_col_const_Multiarray_T get_col_const_Multiarray_c ///< Standard.
#define get_col_Multiarray_T       get_col_Multiarray_c       ///< Standard.

#endif // DPG__templates_multiarray_c_h__INCLUDED
