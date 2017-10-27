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

#ifndef DPG__compute_error_advection_h__INCLUDED
#define DPG__compute_error_advection_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output relating to the linear advection
 *         variables.
 */

struct Simulation;

/** \brief Version of \ref constructor_Error_CE_fptr checking the error of all supported linear advection variables.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_advection_all
	(const struct Simulation* sim ///< Defined for \ref constructor_Error_CE_fptr.
	);

#endif // DPG__compute_error_advection_h__INCLUDED
