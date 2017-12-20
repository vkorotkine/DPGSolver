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

#ifndef DPG__definitions_tol_h__INCLUDED
#define DPG__definitions_tol_h__INCLUDED
/**	\file
 *	\brief Provides the definitions relating to tolerances.
 */

///\{ \name Tolerance values used throughout.
#define EPS      1.0e-15
#define SQRT_EPS 3.162277660168379e-08
#define CX_STEP  1e-30                 ///< Magnitude of the complex step to take for Jacobian evaluation.
///\}

#endif // DPG__definitions_tol_h__INCLUDED
