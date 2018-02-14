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
 *  \brief Provides the interface to templated functions used for parametric geometry processing.
 */

struct Solver_Volume_T;
struct Simulation;

/** \brief Version of \ref constructor_xyz_fptr_T for the parametric cylinder.
 *  \return See brief.
 *
 *  Uses \f$ r-\theta \f$ parametrization to transform from square to circular sections.
 */
const struct const_Multiarray_R* constructor_xyz_cylinder_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_R* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for the n-cube perturbed with trigonometric functions.
 *  \return See brief.
 */
const struct const_Multiarray_R* constructor_xyz_trigonometric_cube_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_R* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);

/** \brief Version of \ref constructor_xyz_fptr_T for a parametric Joukowski to circle blended domain.
 *  \return See brief.
 *
 *  The procedure for the transformation is as follows:
 *  - Find the point corresponding to each input coordinate on the Joukowski and the external cylinder surfaces based on
 *    the reference \f$ r \f$ coordinate.
 *  - Find the final physical coordinate as a linearly blended contribution of the two surface coordinates using
 *    \f$ s \f$ as the blending parameter.
 */
const struct const_Multiarray_R* constructor_xyz_joukowski_parametric_T
	(const char n_type,                      ///< See brief.
	 const struct const_Multiarray_R* xyz_i, ///< See brief.
	 const struct Solver_Volume_T* s_vol,    ///< See brief.
	 const struct Simulation* sim            ///< See brief.
	);
