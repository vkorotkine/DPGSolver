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

#ifndef DPG__definitions_visualization_h__INCLUDED
#define DPG__definitions_visualization_h__INCLUDED
/** \file
 *  \brief Provides the definitions related to the available visualizations.
 */

///\{ \name Available visualization output types.
#define VIS_GEOM_VOLUMES 101 ///< The high-order volume geometry of the \ref Solver_Volume computational elements.
#define VIS_GEOM_EDGES   102 ///< The high-order edge geometry of the \ref Solver_Volume computational elements.
#define VIS_NORMALS      103 ///< The high-order normals of the \ref Solver_Face computational elements.
#define VIS_SOLUTION     104 ///< The high-order solution in the \ref Solver_Volume\*s and \ref Solver_Face\*s.
///\}

#endif // DPG__definitions_visualization_h__INCLUDED
