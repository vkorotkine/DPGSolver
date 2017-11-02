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
/** \file
 */

#include "test_support_math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void z_zyxp_dcc (const int n, const double* x, const double complex* y, double complex* z)
{
	for (int i = 0; i < n; ++i)
		z[i] += y[i]*x[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
