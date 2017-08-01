// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_c_h__INCLUDED
#define DPG__solver_c_h__INCLUDED
/**	\file
 *	\brief Provide solver related functions for complex datatypes.
 *
 *	\todo Potential separate these structs into different files.
 */

#include "containers.h"
#include "containers_c.h"
#include "solver.h"

/// \brief Holds information needed by volumes during the solve.
struct Volume_Solver_c {
	const struct const_Multiarray_d*const xyz,       ///< Moved from S_VOLUME.
	                               *const det_jv_vi, ///< Moved from S_VOLUME.
	                               *const c_vi;      ///< Moved from S_VOLUME.

	struct Multiarray_c*const w_hat,   ///< Solution coefficients.
	                   *const q_hat_v, ///< Local weak gradient coefficients.
	                   *const q_hat,   ///< Weak gradient coefficients.
	                   *const rhs;     ///< Residual.

	struct Volume_Solver_c* next; ///< Pointer to next volume.
};

/// \brief Holds information related to the complex solver context.
struct Context_Solver_c {
	const unsigned int d,     ///< Moved from S_DB.
	                   n_var; ///< Moved from S_DB.

	struct Volume_Solver_c*const volume_head; ///< Pointer to first volume
};

/// \brief Constructor.
struct Context_Solver_c constructor_Context_Solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME_head);

/// \brief Destructor.
void destructor_Context_Solver_c (struct Context_Solver_c* context);

/// \brief Calls function which evaluate q_hat.
extern void compute_GradW_c (const struct S_solver_info*const solver_info, const char stage);

#endif // DPG__solver_c_h__INCLUDED
