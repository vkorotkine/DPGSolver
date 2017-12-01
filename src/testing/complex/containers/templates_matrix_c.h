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

#ifndef DPG__templates_matrix_c_h__INCLUDED
#define DPG__templates_matrix_c_h__INCLUDED
/** \file
 *  \brief Provides the macro definitions used for c-style templating related to the `double` matrix
 *         containers/functions.
 */

#include "templates_matrix_constructors_c.h"
#include "templates_matrix_math_c.h"

///\{ \name Data types
#define Matrix_T       Matrix_c
#define Matrix_R       Matrix_d       ///< 'R'eal.
#define const_Matrix_T const_Matrix_c
#define const_Matrix_R const_Matrix_d ///< 'R'eal.
///\}

///\{ \name Function names
#define set_block_Matrix_T   set_block_Matrix_c
#define set_block_Matrix_T_R set_block_Matrix_c_d
///\}

#endif // DPG__templates_matrix_c_h__INCLUDED
