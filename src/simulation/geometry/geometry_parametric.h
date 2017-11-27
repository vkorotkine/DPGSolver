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

#ifndef DPG__geometry_parametric_h__INCLUDED
#define DPG__geometry_parametric_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for parametric geometry processing.
 */

struct Solver_Volume;
struct Simulation;

/** \brief Version of \ref constructor_xyz_fptr for the parametric cylinder.
 *  \return See brief.
 *
 *  Uses \f$ r-\theta \f$ parametrization to transform from square to circular sections.
 */
const struct const_Multiarray_d* constructor_xyz_cylinder_parametric
	(const struct const_Multiarray_d* xyz_i, ///< See brief.
	 const struct Solver_Volume* s_vol,      ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

#endif // DPG__geometry_parametric_h__INCLUDED
