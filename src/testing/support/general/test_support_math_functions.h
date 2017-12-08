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

#ifndef DPG__test_support_math_functions_h__INCLUDED
#define DPG__test_support_math_functions_h__INCLUDED
/** \file
 *  \brief Provides several standard math functions.
 */

#include "petscmat.h"

/** \brief Computes the relative norm of the difference between the input `Mat` using the infinity norm.
 *  \return See brief. */
double norm_diff_petsc_Mat
	(Mat A0, ///< The Mat for input 0.
	 Mat A1  ///< The Mat for input 1.
	);

#endif // DPG__test_support_math_functions_h__INCLUDED
