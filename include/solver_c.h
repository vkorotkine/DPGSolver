// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_c_h__INCLUDED
#define DPG__solver_c_h__INCLUDED
/**	\file
 *	\brief Provide solver related functions for complex datatypes.
 */

#include "Multiarray.h"
#include "containers_c.h"
#include "Simulation.h"
#include "solver.h"

struct Element_Solver_c;
struct Volume_Solver_c;

/// \brief Struct holding data related to the solve using complex variables.
struct Solver_c {
	const unsigned int d,     ///< Defined in \ref Simulation.
	                   n_var; ///< Defined in \ref Simulation.

	struct Element_Solver_c*const element_head; ///< Pointer to the head of the \ref Element_Solver_c list.
	struct Volume_Solver_c*       volume_head;  ///< Pointer to the head of the \ref Volume_Solver_c list.
};

/// \brief Struct holding data related to the Elements for the solve using complex variables.
struct Element_Solver_c {

	struct Element_Solver_c* prev, ///< Pointer to the previous \ref Element_Solver_c in the list.
	                       * next; ///< Pointer to the next     \ref Element_Solver_c in the list.
};

/// \brief Struct holding data related to the Volumes for the solve using complex variables.
struct Volume_Solver_c {
	const struct const_Multiarray_d*const xyz,       ///< Moved from \ref S_VOLUME.
	                               *const det_jv_vi, ///< Moved from \ref S_VOLUME.
	                               *const c_vi;      ///< Moved from \ref S_VOLUME.

	struct Multiarray_c*const w_hat,   ///< Solution coefficients.
	                   *const q_hat_v, ///< Volume contribution to the weak gradient coefficients.
	                   *const q_hat,   ///< Weak gradient coefficients.
	                   *const rhs;     ///< Residual.

	const struct Element_Solver_c*const element; ///< Pointer to the element associated with this volume.

	struct Volume_Solver_c* prev, ///< Pointer to previous \ref Volume_Solver_c in the list.
	                      * next; ///< Pointer to next     \ref Volume_Solver_c in the list.
};

/// \brief Constructor for \ref Solver_c.
struct Solver_c* constructor_Solver_c
	(const struct Simulation*const simulation ///< Standard.
	);

/// \brief Destructor for \ref Solver_c.
void destructor_Solver_c
	(struct Solver_c* solver_c ///< Standard.
	);




/// \brief Calls function which evaluate q_hat.
/// \todo update this.
extern void compute_GradW_c (const struct S_solver_info*const solver_info, const char stage);

#endif // DPG__solver_c_h__INCLUDED
