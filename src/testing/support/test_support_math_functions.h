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

#include <complex.h>

/// \brief Variation on the axpy (BLAS 1) function: z = y*x + z (`double`,`double complex`,`double complex`).
void z_yxpz_dcc
	(const int n,     ///< The number of entries.
	 const double* x, ///< Input x.
	 const double complex* y, ///< Input y.
	 double complex* z        ///< Location to store the sum.
	);

#endif // DPG__test_support_math_functions_h__INCLUDED
