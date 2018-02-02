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

#ifndef DPG__compute_error_navier_stokes_h__INCLUDED
#define DPG__compute_error_navier_stokes_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output relating to the Euler variables.
 */

struct Simulation;

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of velocity and temperature.
 *  \return See brief.
 *
 *  This function should be used when the exact solution is available.
 */
struct Error_CE* constructor_Error_CE_navier_stokes_uvwt
	(const struct Simulation* sim ///< Defined for \ref constructor_Error_CE_fptr.
	);

#endif // DPG__compute_error_navier_stokes_h__INCLUDED
