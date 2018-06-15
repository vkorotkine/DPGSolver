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
 *  \brief Provides templated containers and functions relating to boundary conditions for the Euler equation.
 */

#include "def_templates_boundary.h"
#include "def_templates_face_solver.h"

struct Boundary_Value_Input_T;
struct Boundary_Value_T;
struct Solver_Face_T;
struct Simulation;

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the Riemann invariant values.
 *
 *  Reference: (Section 2.2, \cite Carlson2011); note the typo in eq. (14).
 */
void constructor_Boundary_Value_T_euler_riemann
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the slip wall values.
 *
 *  The slip wall boundary condition sets the ghost state variables as follows:
 *  - density:  equal to internal density;
 *  - pressure: equal to internal pressure;
 *  - velocity: equal to internal velocity reflected along the tangent to the boundary.
 *
 *  \note As the density, pressure and velocity magnitude of the ghost state are all equal to those of the internal
 *        state, it is immediately noted that the total energy of the ghost state is also equal to the internal value.
 */
void constructor_Boundary_Value_T_euler_slipwall
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the external upwind values.
void constructor_Boundary_Value_T_euler_supersonic_inflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/// \brief Version of \ref constructor_Boundary_Value_fptr_T computing members using the internal interpolated values.
void constructor_Boundary_Value_T_euler_supersonic_outflow
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members by imposing the back pressure.
 *
 *  Reference: (Section 2.4, \cite Carlson2011).
 */
void constructor_Boundary_Value_T_euler_back_pressure
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

/** \brief Version of \ref constructor_Boundary_Value_fptr_T computing members by imposing the total temperature and
 *         pressure.
 *
 *  Eqs. (38/47) (\cite Carlson2011) imply that the velocity should be normal to the boundary.
 *
 *  Reference: (Section 2.7, \cite Carlson2011); (Eqs. (3.9), (8.58), \cite Toro2009).
 */
void constructor_Boundary_Value_T_euler_total_tp
	(struct Boundary_Value_T* bv,               ///< See brief.
	 const struct Boundary_Value_Input_T* bv_i, ///< See brief.
	 const struct Solver_Face_T* face,          ///< See brief.
	 const struct Simulation* sim               ///< See brief.
	);

#include "undef_templates_boundary.h"
#include "undef_templates_face_solver.h"
