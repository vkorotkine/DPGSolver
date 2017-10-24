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

#ifndef DPG__solution_h__INCLUDED
#define DPG__solution_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for solution specification (initialization).
 */

struct Simulation;
struct Solver_Volume;
struct Solver_Face;

/** \brief Function pointer to volume solution setting function.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Solver_Volume.
 */
typedef void (*set_sol_coef_v_fptr)
	(const struct Simulation* sim,
	 struct Solver_Volume* volume
	);

/** \brief Function pointer to face solution setting function.
 *  \param sim  \ref Simulation.
 *  \param face \ref Solver_Face.
 */
typedef void (*set_sol_coef_f_fptr)
	(const struct Simulation* sim,
	 struct Solver_Face* face
	);

/** \brief Function pointer to the function setting the source contribution of \ref Solver_Volume::rhs.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Solver_Volume.
 */
typedef void (*compute_source_fptr)
	(const struct Simulation* sim,
	 struct Solver_Volume* volume
	);

// Interface functions ********************************************************************************************** //

/** \brief Set up the initial solution for the simulation. Sets:
 *	- \ref Solver_Volume::sol_coef;
 *	- \ref Solver_Volume::grad_coef (if applicable);
 *	- \ref Solver_Face::sol_coef    (if applicable);
 *	- \ref Solver_Face::grad_coef   (if applicable).
 */
void set_initial_solution
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Function pointer to be used for \ref Test_Case::set_grad_coef_v when there are no volume solution gradients.
void set_grad_coef_v_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref set_sol_coef_v_fptr.
	 struct Solver_Volume* volume  ///< Defined for \ref set_sol_coef_v_fptr.
	);

/// \brief Function pointer to be used for \ref Test_Case::set_sol_coef_f when there is no face solution.
void set_sol_coef_f_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref set_sol_coef_f_fptr.
	 struct Solver_Face* face      ///< Defined for \ref set_sol_coef_f_fptr.
	);

/// \brief Function pointer to be used for \ref Test_Case::set_grad_coef_f when there are no face solution gradients.
void set_grad_coef_f_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref set_sol_coef_f_fptr.
	 struct Solver_Face* face      ///< Defined for \ref set_sol_coef_f_fptr.
	);

/** \brief Contructor for a \ref const_Multiarray_d\* holding the xyz coordinates at the volume solution nodes.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_xyz_vs
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume* volume  ///< \ref Solver_Volume.
	);

/// \brief Function pointer to be used for \ref Test_Case::compute_source when there is no source term.
void compute_source_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref compute_source_fptr.
	 struct Solver_Volume* volume  ///< Defined for \ref compute_source_fptr.
	);

#endif // DPG__solution_h__INCLUDED
