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
 *  \brief Provides templated functions used for computing the volume contributions to the right and left-hand side
 *         (rlhs) terms of supported schemes.
 */

#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_volume_solver.h"
#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_flux.h"

struct Flux_Input_T;
struct Flux_Ref_T;
struct Simulation;
struct Solver_Storage_Implicit;
struct Volume;
struct Solver_Volume_T;

/** \brief Function pointer to the constructor for the solution evaluated at the volume cubature nodes.
 *  \return Standard.
 *
 *  \param s_vol The current volume.
 */
typedef const struct const_Multiarray_T* (*constructor_sol_vc_fptr_T)
	(const struct Solver_Volume_T*const s_vol
	);

/** \brief Function pointer to the destructor for the solution evaluated at the volume cubature nodes.
 *
 *  \param sol_vc To be destructed.
 */
typedef void (*destructor_sol_vc_fptr_T)
	(const struct const_Multiarray_T* sol_vc
	);

/** \brief Pointer to the function used to evaluate the rhs (and optionally lhs) volume terms.
 *
 *  \param flux_r \ref Flux_Ref_T.
 *  \param s_vol  \ref Solver_Volume_T.
 *  \param ssi    \ref Solver_Storage_Implicit.
 */
typedef void (*compute_rlhs_v_fptr_T)
	(const struct Flux_Ref_T*const flux_r,
	 struct Solver_Volume_T*const s_vol,
	 struct Solver_Storage_Implicit*const ssi
	);

/// \brief Container for volume solver parameters related to con/de'structor's.
struct S_Params_Volume_Structor_T {
	constructor_sol_vc_fptr_T constructor_sol_vc;  ///< Pointer to the appropriate function.
	constructor_sol_vc_fptr_T constructor_grad_vc; ///< Pointer to the appropriate function.

	destructor_sol_vc_fptr_T  destructor_sol_vc;   ///< Pointer to the appropriate function.
	destructor_sol_vc_fptr_T  destructor_grad_vc;  ///< Pointer to the appropriate function.
};

/// \brief Container for the reference flux related parameters.
struct Flux_Ref_T {
	const struct const_Multiarray_T* fr;       ///< The reference flux.
	const struct const_Multiarray_T* dfr_ds;   ///< The reference flux Jacobians with respect to the solution.
	const struct const_Multiarray_T* dfr_dg;   ///< The reference flux Jacobians with respect to the solution gradients.
	const struct const_Multiarray_T* d2fr_ds2; ///< The reference flux Hessians  with respect to the solution.
};

// Interface functions ********************************************************************************************** //

/// \brief Set the parameters of \ref S_Params_Volume_Structor_T.
void set_S_Params_Volume_Structor_T
	(struct S_Params_Volume_Structor_T* spvs, ///< \ref S_Params_Volume_Structor_T.
	 const struct Simulation* sim             ///< \ref Simulation.
	);

/** \brief Call \ref constructor_Flux_Ref using data from the current volume.
 *  \return See brief. */
struct Flux_Ref_T* constructor_Flux_Ref_vol_T
	(const struct S_Params_Volume_Structor_T* spvs, ///< \ref S_Params_Volume_Structor_T.
	 struct Flux_Input_T* flux_i,                   ///< \ref Flux_Input_T.
	 const struct Solver_Volume_T* s_vol,           ///< \ref Solver_Volume_T.
	 const struct Simulation* sim                   ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Flux_Ref_T container.
void destructor_Flux_Ref_T
	(struct Flux_Ref_T* flux_ref ///< Standard.
	);

/** \brief Version of \ref compute_rlhs_v_fptr_T computing only the rhs term.
 *
 *  The function is "dg-like" in the sense that it computes the rhs term as for the DG scheme, but can be used with test
 *  space not equal to trial space (Petrov-Galerkin).
 */
void compute_rhs_v_dg_like_T
	(const struct Flux_Ref_T*const flux_r,    ///< See brief.
	 struct Solver_Volume_T*const s_vol,      ///< See brief.
	 struct Solver_Storage_Implicit*const ssi ///< See brief.
	);

/** \brief Constructor for the lhs volume term of 1st order equations only (i.e. flux having dependence only on \ref
 *         Solver_Volume_T::sol_coef.
 *  \return See brief. */
struct Matrix_T* constructor_lhs_v_1_T
	(const struct Flux_Ref_T*const flux_r,    ///< \ref Flux_Ref_T.
	 const struct Solver_Volume_T*const s_vol ///< \ref Solver_Volume_T.
	);

/** \brief Constructor for the partial lhs volume term of 2nd order equations only (i.e. flux having dependence only on
 *         \ref Solver_Volume_T::grad_coef).
 *  \return See brief.
 *
 *  The lhs contribution computed here is termed 'p'artial as it must subsequently be multiplied by the
 *  d_g_coef__d_s_coef terms in order to obtain the final lhs terms.
 */
struct Matrix_T* constructor_lhs_p_v_2_T
	(const struct Flux_Ref_T*const flux_r,    ///< \ref Flux_Ref_T.
	 const struct Solver_Volume_T*const s_vol ///< \ref Solver_Volume_T.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vr_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vr_vc_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vs_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vs_vc_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vt_vc operator.
 *  \return See brief. */
const struct Operator* get_operator__cv0_vt_vc_T
	(const struct Solver_Volume_T*const s_vol ///< The current volume.
	);

/** \brief Get the appropriate sub-range of the \ref Solver_Element::tw1_vt_vc operators.
 *  \return See brief. */
struct Multiarray_Operator get_operator__tw1_vt_vc_T
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Get the appropriate sub-range of the \ref Solver_Element::cv1_vt_vc operators.
 *  \return See brief. */
struct Multiarray_Operator get_operator__cv1_vt_vc_T
	(const struct Solver_Volume_T*const s_vol ///< The current volume.
	);

#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_flux.h"
