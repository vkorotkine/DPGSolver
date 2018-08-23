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

#ifndef DPG__definitions_numerical_flux_h__INCLUDED
#define DPG__definitions_numerical_flux_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to numerical fluxes.
 */

///\{ \name Definitions for the available numerical fluxes.
#define NUM_FLUX_INVALID   0

#define NUM_FLUX_UPWIND    11

#define NUM_FLUX_BR2_STABLE 21
#define NUM_FLUX_CDG2       22

#define NUM_FLUX_ROE_PIKE       31
#define NUM_FLUX_LAX_FRIEDRICHS 32
///\}

#endif // DPG__definitions_numerical_flux_h__INCLUDED
