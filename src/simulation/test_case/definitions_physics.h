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

#ifndef DPG__definitions_physics_h__INCLUDED
#define DPG__definitions_physics_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to physics.
 */

///\{ \name Definitions for ratio of specific heats.
#define GAMMA  1.4
#define GM1    0.4
#define GM3   -1.6
///\}

///\{ \name Minimum permited value for physical quantities which should be positive.
#define EPS_PHYS 1.0e-13
///\}

///\{ \name Definitions related to the method in which the viscosity should be computed.
#define VISCOSITY_INVALID    -1  ///< Invalid value.
#define VISCOSITY_CONSTANT   101 ///< Constant viscosity.
#define VISCOSITY_SUTHERLAND 102 ///< Viscosity computed using Sutherland's formula.
///\}

#endif // DPG__definitions_physics_h__INCLUDED
