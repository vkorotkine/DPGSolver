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

#ifndef DPG__definitions_solution_h__INCLUDED
#define DPG__definitions_solution_h__INCLUDED
/** \file
 *  \brief Provides the definitions relating to solutions.
 */

///\{ \name Definitions for the available options for advection velocity computing functions.
#define ADVECTION_TYPE_CONST       1 ///< See \ref compute_b_adv_constant.
#define ADVECTION_TYPE_VORTEX      2 ///< See \ref compute_b_adv_vortex.
///\}

///\{ \name Definitions for the available options for boundary solution perturbation.
#define BOUNDARY_PERTURB_TYPE_NONE   1 ///< No perturbation.
#define BOUNDARY_PERTURB_TYPE_TRIG_X 2 ///< Perturb with trigonometric function on x-coordinate boundaries.
///\}

#endif // DPG__definitions_solution_h__INCLUDED
