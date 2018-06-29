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

#ifndef DPG__optimize_h__INCLUDED
#define DPG__optimize_h__INCLUDED

struct Simulation;
struct Optimization_Case;


/** \brief Optimize the the geometry to minimize a specified functional.
 */
void optimize(
	struct Simulation* sim ///< Standard (the simulation data structure)
	);


void output_NURBS_patch_information(struct Optimization_Case* optimization_case);

#endif // DPG__optimize_h__INCLUDED