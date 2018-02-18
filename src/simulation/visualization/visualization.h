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

#ifndef DPG__visualization_h__INCLUDED
#define DPG__visualization_h__INCLUDED
/** \file
 *  \brief Provides functions used for visualization of outputs.
 *
 *  In preparation for parallelization of the code, a 'p'arallel and 's'erial file is generated for each of the outputs.
 *  For example, for the unstructured Paraview format, both '.pvtu' and '.vtu' files are generated.
 */

struct Simulation;

/// \brief Output the visualization of the specified output.
void output_visualization
	(struct Simulation* sim, ///< \ref Simulation.
	 const int vis_type      ///< The type of visualization. Options: see \ref definitions_visualization.h.
	);

#endif // DPG__visualization_h__INCLUDED
