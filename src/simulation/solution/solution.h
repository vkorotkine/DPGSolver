// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solution_h__INCLUDED
#define DPG__solution_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for solution specification.
 */

struct Simulation;
struct Intrusive_List;

/** \brief Set up the initial solution for the simulation. Computes:
 *	- \ref Volume::sol_coef;
 *	- \ref Volume::grad_coef (if applicable).
 */
void set_up_solution
	(struct Simulation* sim,               ///< \ref Simulation.
	 struct Intrusive_List* solver_volumes ///< The \ref Solver_Volume list for which to set up the solution.
	);

#endif // DPG__solution_h__INCLUDED
