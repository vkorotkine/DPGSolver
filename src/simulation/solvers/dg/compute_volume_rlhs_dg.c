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
 */

#include "compute_volume_rlhs_dg.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"

#include "volume_solver_dg.h"
#include "element_solver_dg.h"

#include "intrusive.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the solution evaluated at the volume cubature nodes.
 *  \return Standard.
 *
 *  \param volume The current volume.
 *  \param sim    \ref Simulation.
 */
typedef const struct const_Multiarray_d* (*constructor_sol_vc_fptr)
	(struct Volume* volume,
	 const struct Simulation* sim
	);

/** \brief Destructor for the solution evaluated at the volume cubature nodes.
 *
 *  \param sol_vc To be destructed.
 */
typedef void (*destructor_sol_vc_fptr)
	(const struct const_Multiarray_d* sol_vc
	);

/// \brief Container for function pointers.
struct F_Ptrs {
	constructor_sol_vc_fptr constructor_sol_vc; ///< Pointer to the appropriate function.
	destructor_sol_vc_fptr  destructor_sol_vc;  ///< Pointer to the appropriate function.
};

/** \brief Set the function pointers in \ref F_Ptrs.
 *  \return A statically allocated \ref F_Ptrs container. */
struct F_Ptrs set_f_ptrs
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set the memory of the rhs and lhs (if applicable) terms to zero for the volumes.
static void zero_memory_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rhs_dg (const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);

	struct F_Ptrs f_ptrs = set_f_ptrs(sim);

/// \todo Check if this is required.
	zero_memory_volumes(sim);

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume*        vol   = (struct Volume*) curr;
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		const struct const_Multiarray_d* sol_vc = f_ptrs.constructor_sol_vc(vol,sim);

print_Multiarray_d(s_vol->sol_coef);
print_const_Multiarray_d(sol_vc);

		const ptrdiff_t n_vc = sol_vc->extents[0];
		const ptrdiff_t extents[3] = { n_vc,d,n_eq };

// Make a constructor_flux function:
//    - Should set function pointer to the appropriate function in the test_case construction.
//    - Since they are already split by dimension, make each of these functions separate and use a compile-time constant
//    for the number of variables, dimension.
//    - Should take inputs of: sol, grad_sol, sim.
//    - Ensure that both fluxes use '+=' instead of the '=' they are currently using after zeroing initial memory.
//    - Should return fr_vc, dfrds_vc, ... (whatever is needed based on test_case parameters)
//    - Should sum the inviscid and viscous contributions.
		struct Multiarray_d* f_vc = constructor_zero_Multiarray_d('C',3,extents);
EXIT_UNSUPPORTED;

		f_ptrs.destructor_sol_vc(sol_vc);

		destructor_const_Multiarray_d(f_vc);
	}

	EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the solution evaluated at the volume cubature nodes using interpolation.
 *  \return Standard. */
static const struct const_Multiarray_d* constructor_sol_vc_interp
	(struct Volume* volume,       ///< The current volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for the solution evaluated at the volume cubature nodes using collocation.
 *  \return Standard. */
static const struct const_Multiarray_d* constructor_sol_vc_col
	(struct Volume* volume,       ///< The current volume.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for the return value of \ref constructor_sol_vc_interp.
static void destructor_sol_vc_interp
	(const struct const_Multiarray_d* sol_vc ///< To be destructed.
	);

/// \brief Destructor for the return value of \ref constructor_sol_vc_col.
static void destructor_sol_vc_col
	(const struct const_Multiarray_d* sol_vc ///< To be destructed.
	);

struct F_Ptrs set_f_ptrs (const struct Simulation* sim)
{
	struct F_Ptrs f_ptrs;

	if (!sim->collocated) {
		f_ptrs.constructor_sol_vc = constructor_sol_vc_interp;
		f_ptrs.destructor_sol_vc  = destructor_sol_vc_interp;
	} else {
		f_ptrs.constructor_sol_vc = constructor_sol_vc_col;
		f_ptrs.destructor_sol_vc  = destructor_sol_vc_col;
	}

	return f_ptrs;
}

static void zero_memory_volumes (const struct Simulation* sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		set_to_value_Multiarray_d(((struct DG_Solver_Volume*)curr)->rhs,0.0);
	}
}

// Level 0 ********************************************************************************************************** //

static const struct const_Multiarray_d* constructor_sol_vc_interp (struct Volume* volume, const struct Simulation* sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Volume* s_volume = (struct Solver_Volume*) volume;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) volume->element;

	const int p = s_volume->p_ref;
	const struct Operator* cv0_vs_vc =
		(!volume->curved ? get_Multiarray_Operator(e->cv0_vs_vcs,(ptrdiff_t[]){0,0,p,p})
		                 : get_Multiarray_Operator(e->cv0_vs_vcc,(ptrdiff_t[]){0,0,p,p}) );

	const struct const_Multiarray_d* s_coef = (const struct const_Multiarray_d*) s_volume->sol_coef;

	return constructor_mm_NN1_Operator_const_Multiarray_d(cv0_vs_vc,s_coef,'C',op_format,s_coef->order,NULL);
}

static const struct const_Multiarray_d* constructor_sol_vc_col (struct Volume* volume, const struct Simulation* sim)
{
	UNUSED(sim);

	struct Solver_Volume* s_volume = (struct Solver_Volume*) volume;
	return (const struct const_Multiarray_d*) s_volume->sol_coef;
}

static void destructor_sol_vc_interp (const struct const_Multiarray_d* sol_vc)
{
	destructor_const_Multiarray_d(sol_vc);
}

static void destructor_sol_vc_col (const struct const_Multiarray_d* sol_vc)
{
	UNUSED(sol_vc);
}
