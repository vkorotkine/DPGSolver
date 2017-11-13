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

#ifndef DPG__compute_volume_rlhs_h__INCLUDED
#define DPG__compute_volume_rlhs_h__INCLUDED
/** \file
 *  \brief Provides functions used for computing the volume contributions to the right and left-hand side (rlhs) terms
 *         of supported schemes.
 */

struct Flux_Input;
struct Simulation;
struct Volume;
struct Solver_Volume;

/** \brief Function pointer to the constructor for the solution evaluated at the volume cubature nodes.
 *  \return Standard.
 *
 *  \param s_vol The current volume.
 *  \param sim   \ref Simulation.
 */
typedef const struct const_Multiarray_d* (*constructor_sol_vc_fptr)
	(const struct Solver_Volume* s_vol,
	 const struct Simulation* sim
	);

/** \brief Function pointer to the destructor for the solution evaluated at the volume cubature nodes.
 *
 *  \param sol_vc To be destructed.
 */
typedef void (*destructor_sol_vc_fptr)
	(const struct const_Multiarray_d* sol_vc
	);

/// \brief Container for volume solver parameters related to con/de'structor's.
struct S_Params_Volume_Structor {
	constructor_sol_vc_fptr constructor_sol_vc; ///< Pointer to the appropriate function.
	destructor_sol_vc_fptr  destructor_sol_vc;  ///< Pointer to the appropriate function.
};

/// \brief Container for the reference flux related parameters.
struct Flux_Ref {
	const struct const_Multiarray_d* fr;     ///< The reference flux.
	const struct const_Multiarray_d* dfr_ds; ///< The reference flux Jacobians with respect to the solution.
	const struct const_Multiarray_d* dfr_dg; ///< The reference flux Jacobians with respect to the solution gradients.
};

// Interface functions ********************************************************************************************** //

/** \brief Set the parameters of \ref S_Params_Volume_Structor.
 *  \return A statically allocated \ref S_Params container. */
void set_S_Params_Volume_Structor
	(struct S_Params_Volume_Structor* spvs, ///< \ref S_Params_Volume_Structor.
	 const struct Simulation* sim           ///< \ref Simulation.
	);

/** \brief Call \ref constructor_Flux_Ref using data from the current volume.
 *  \return See brief. */
struct Flux_Ref* constructor_Flux_Ref_vol
	(const struct S_Params_Volume_Structor* spvs, ///< \ref S_Params_Volume_Structor.
	 struct Flux_Input* flux_i,                   ///< \ref Flux_Input.
	 const struct Solver_Volume* s_vol,           ///< \ref Solver_Volume.
	 const struct Simulation* sim                 ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Flux_Ref container.
void destructor_Flux_Ref
	(struct Flux_Ref* flux_ref ///< Standard.
	);

/** \brief Constructor for the lhs volume term of 1st order equations only.
 *  \return See brief. */
struct Matrix_d* constructor_lhs_v_1
	(const struct Flux_Ref* flux_r,     ///< \ref Flux_Ref.
	 const struct Solver_Volume* s_vol, ///< \ref Solver_Volume.
	 const struct Simulation* sim       ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vs_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vs_vc
	(const struct Solver_Volume* s_vol ///< The current volume.
	);

/** \brief Get the appropriate sub-range of the \ref Solver_Element::tw1_vt_vc operators.
 *  \return See brief. */
struct Multiarray_Operator get_operator__tw1_vt_vc__rlhs
	(const struct Solver_Volume* s_vol ///< The current volume.
	);

#endif // DPG__compute_volume_rlhs_h__INCLUDED
