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
/** \file
 *  \brief Provides the definitions related to the boundary conditions.
 *
 *  The naming convention for the periodic boundary conditions omits the 'BC_' prefix as they are only used for
 *  preprocessing and not as actual boundary conditions.
 *
 *  The "ALT"ernate values for the same boundary conditions are required such that volumes having multiple vertices on a
 *  boundary are not identified as boundary volumes when these vertices are not on the same boundary. The alternate
 *  values are thus used to distinguish between separate portions of the boundary employing the same boundary condition.
 *
 *  \note The value of these constants must be identical to those in the $ROOT/input/meshes/parameters.geo file.
 */

#define BC_INVALID -1 ///< Value to indicate an invalid boundary condition.

#define BC_STEP_SC      10000 ///< The step in the boundary conditions between straight and curved BCs.
#define BC_CURVED_START 20000 ///< The value after which boundaries are considered to be curved.

///\{ \name The values for the various boundary conditions after taking the modulus with \ref BC_STEP_SC.

// Advection
#define BC_INFLOW       1
#define BC_INFLOW_ALT1  2
#define BC_INFLOW_ALT2  3
#define BC_OUTFLOW      11
#define BC_OUTFLOW_ALT1 12
#define BC_OUTFLOW_ALT2 13

// Diffusion
#define BC_DIRICHLET      21
#define BC_DIRICHLET_ALT1 22
#define BC_NEUMANN        31
#define BC_NEUMANN_ALT1   32

// Euler
#define BC_RIEMANN        101
#define BC_SLIPWALL       102
#define BC_BACKPRESSURE   103
#define BC_TOTAL_TP       104
#define BC_SUPERSONIC_IN  105
#define BC_SUPERSONIC_OUT 106

// Navier-Stokes
#define BC_NOSLIP_ADIABATIC    111
#define BC_NOSLIP_DIABATIC     112
#define BC_NOSLIP_ALL_ROTATING 113

///\}

///\{ \name Options for specialized viscous boundary conditions.
#define VISCOUS_BC_INVALID 101 ///< Invalid value for this set of options.

#define DIABATIC_FLUX_CONSTANT_ZERO 201 ///< Constant value (equal to 0.0) for the viscous energy flux.
#define DIABATIC_FLUX_CONSTANT      202 ///< Constant value for the viscous energy flux.

#define NO_SLIP_ROTATING 301 ///< Obtain the no-slip conditions using the angular velocity of a rotating cylinder.

/// Obtain all of the remaining BCs (i.e. all excluding velocity) from the density (\f$ \rho \f$) and total energy (E).
#define NO_SLIP_ALL__RHO_E 401

#define NO_SLIP_ALL_ROTATING_RHO_E 1001 ///< Compound \ref NO_SLIP_ROTATING and \ref NO_SLIP_ALL__RHO_E.
///\}

#endif // DPG__definitions_bc_h__INCLUDED
