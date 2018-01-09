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
 *  \brief Provides the interface to functions used for solution specification (initialization).
 */

struct const_Multiarray_R;
struct const_Multiarray_T;
struct Simulation;
struct Solver_Volume_T;
struct Solver_Face_T;
struct Solution_Container_T;

/** \brief Function pointer to a function which constructs the solution given the input xyz coordinates.
 *  \return A \ref const_Multiarray_T\* container holding the solution.
 *
 *  \param xyz Input xyz coordinates.
 *  \param sim \ref Simulation.
 */
typedef const struct const_Multiarray_T* (*constructor_sol_fptr_T)
	(const struct const_Multiarray_R* xyz,
	 const struct Simulation* sim
	);

/** \brief `mutable` version of \ref constructor_sol_fptr_T.
 *  \return See brief.
 *
 *  \param xyz See brief.
 *  \param sim See brief.
 */
typedef struct Multiarray_T* (*mutable_constructor_sol_fptr_T)
	(const struct const_Multiarray_R* xyz,
	 const struct Simulation* sim
	);

/** \brief Function pointer to volume solution setting function.
 *  \param sim      \ref Simulation.
 *  \param sol_cont \ref Solution_Container_T.
 */
typedef void (*set_sol_fptr_T)
	(const struct Simulation* sim,
	 struct Solution_Container_T sol_cont
	);

/** \brief Function pointer to the function setting the source contribution of the rhs term.
 *  \param sim   \ref Simulation.
 *  \param s_vol \ref Solver_Volume_T.
 *  \param rhs   Memory in which to add the rhs contribution.
 */
typedef void (*compute_source_rhs_fptr_T)
	(const struct Simulation* sim,
	 const struct Solver_Volume_T* s_vol,
	 struct Multiarray_T* rhs
	);

/// Container for members relating to the solution computation.
struct Solution_Container_T {
	const char ce_type,   ///< The type of computational element associated with the solution data being set.
	           cv_type,   ///< The format in which to return the solution. Options: 'c'oefficients, 'v'alues.
	           node_kind; ///< The kind of nodes to be used. Options: 's'olution, 'c'ubature.

	struct Solver_Volume_T* volume; ///< \ref Solver_Volume_T.
	struct Solver_Face_T* face;     ///< \ref Solver_Face_T.

	struct Multiarray_T* sol; ///< The container for the computed solution.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_sol_fptr_T to be used what a call to this function should be invalid.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_sol_invalid_T
	(const struct const_Multiarray_R* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/** \brief Set up the initial solution for the simulation. Sets:
 *	- \ref Solver_Volume_T::sol_coef;
 *	- \ref Solver_Volume_T::grad_coef (if applicable);
 *	- \ref Solver_Face_T::nf_coef     (if applicable);
 */
void set_initial_solution_T
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Function pointer to be used for \ref Test_Case_T::set_sol or \ref Test_Case_T::set_grad when this solution is not
 *         required. */
void set_sg_do_nothing_T
	(const struct Simulation* sim,      ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief Contructor for a \ref const_Multiarray_T\* holding the xyz coordinates associated with the
 *         \ref Solution_Container_T.
 *  \return See brief. */
const struct const_Multiarray_R* constructor_xyz_sol_T
	(const struct Simulation* sim, ///< \ref Simulation.
	 const struct Solution_Container_T* sol_cont ///< \ref Solution_Container_T.
	);

/** \brief Contructor for a \ref Multiarray_T\* holding the solution at volume nodes of input kind.
 *  \return See brief. */
struct Multiarray_T* constructor_sol_v_T
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume_T* s_vol, ///< \ref Solver_Volume_T.
	 const char node_kind           ///< The kind of node. Options: 's'olution, 'c'ubature.
	);

/// \brief Function pointer to be used for \ref Test_Case_T::compute_source_rhs when there is no source term.
void compute_source_rhs_do_nothing_T
	(const struct Simulation* sim,      ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs           ///< See brief.
	);

/// \brief Update \ref Solution_Container_T::sol based on the input solution values.
void update_Solution_Container_sol_T
	(struct Solution_Container_T*const sol_cont, ///< Defined for \ref set_sol_fptr_T.
	 struct Multiarray_T*const sol,              ///< The solution values.
	 const struct Simulation*const sim           ///< \ref Simulation.
	);

/** \brief Constructor for the xyz coordinates evaluated at the volume cubature nodes using interpolation.
 *  \return See brief. */
const struct const_Multiarray_R* constructor_xyz_vc_interp_T
	(const struct Solver_Volume_T* s_vol, ///< The current volume.
	 const struct Simulation* sim       ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::tw0_vt_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__tw0_vt_vc_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);
