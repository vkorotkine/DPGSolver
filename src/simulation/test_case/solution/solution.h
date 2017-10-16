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
 *  \brief Provides the interface to functions used for solution computation.
 */

struct Simulation;
struct Solver_Volume;
struct Solver_Face;

/** \brief Function pointer to volume solution computation function.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Solver_Volume.
 */
typedef void (*compute_sol_coef_v_fptr)
	(const struct Simulation* sim,
	 struct Solver_Volume* volume
	);

/** \brief Function pointer to face solution computation function.
 *  \param sim  \ref Simulation.
 *  \param face \ref Solver_Face.
 */
typedef void (*compute_sol_coef_f_fptr)
	(const struct Simulation* sim,
	 struct Solver_Face* face
	);

// Interface functions ********************************************************************************************** //

/// \brief Compute the solution for the current simulation.
void compute_solution
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Function pointer to be used for compute_grad_coef_v when there are no volume solution gradients.
void compute_grad_coef_v_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref compute_sol_coef_v_fptr.
	 struct Solver_Volume* volume  ///< Defined for \ref compute_sol_coef_v_fptr.
	);

/// \brief Function pointer to be used for compute_sol_coef_f when there is no face solution.
void compute_sol_coef_f_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref compute_sol_coef_f_fptr.
	 struct Solver_Face* face      ///< Defined for \ref compute_sol_coef_f_fptr.
	);

/// \brief Function pointer to be used for compute_grad_coef_f when there are no face solution gradients.
void compute_grad_coef_f_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref compute_sol_coef_f_fptr.
	 struct Solver_Face* face      ///< Defined for \ref compute_sol_coef_f_fptr.
	);

/** \brief Contructor for a \ref const_Multiarray_d\* holding the xyz coordinates at the volume solution nodes.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_xyz_vs
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume* volume  ///< \ref Solver_Volume.
	);

#endif // DPG__solution_h__INCLUDED
