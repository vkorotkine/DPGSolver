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

#ifndef DPG__definitions_bc_h__INCLUDED
#define DPG__definitions_bc_h__INCLUDED
/**	\file
 *	\brief Provides the definitions related to the boundary conditions.
 *
 *	The naming convention for the periodic boundary conditions omits the 'BC_' prefix as they are only used for
 *	preprocessing and not as actual boundary conditions.
 *
 *	\note The value of these constants must be identical to those in the $ROOT/input/meshes/parameters.geo file.
 */

#define BC_INVALID -1 ///< Value to indicate an invalid boundary condition.

#define BC_STEP_SC      10000 ///< The step in the boundary conditions between straight and curved BCs.
#define BC_CURVED_START 20000 ///< The value after which boundaries are considered to be curved.

///\{ \name The values for the various boundary conditions after taking the modulus with \ref BC_STEP_SC.

// General
#define PERIODIC_XL 51
#define PERIODIC_XR 52
#define PERIODIC_YL 53
#define PERIODIC_YR 54
#define PERIODIC_ZL 55
#define PERIODIC_ZR 56

// Advection
#define BC_INFLOW  1
#define BC_OUTFLOW 2

// Diffusion
#define BC_DIRICHLET 11
#define BC_NEUMANN   12

// Euler
#define BC_RIEMANN        101
#define BC_SLIPWALL       102
#define BC_BACKPRESSURE   103
#define BC_TOTAL_TP       104
#define BC_SUPERSONIC_IN  105
#define BC_SUPERSONIC_OUT 106

// Navier-Stokes
#define BC_NOSLIP_T         111
#define BC_NOSLIP_T_ALL     112
#define BC_NOSLIP_ADIABATIC 116
#define BC_NOSLIP_DIABATIC  117
///\}


#endif // DPG__definitions_bc_h__INCLUDED
