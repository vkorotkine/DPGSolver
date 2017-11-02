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

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "volume_solver_dg.h"
#include "element_solver_dg.h"

#include "flux.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

struct Flux_Ref;

/** \brief Function pointer to the constructor for the solution evaluated at the volume cubature nodes.
 *  \return Standard.
 *
 *  \param volume The current volume.
 *  \param sim    \ref Simulation.
 */
typedef const struct const_Multiarray_d* (*constructor_sol_vc_fptr)
	(struct Volume* volume,
	 const struct Simulation* sim
	);

/** \brief Function pointer to the destructor for the solution evaluated at the volume cubature nodes.
 *
 *  \param sol_vc To be destructed.
 */
typedef void (*destructor_sol_vc_fptr)
	(const struct const_Multiarray_d* sol_vc
	);

/** \brief Function pointer to the function used to evaluate the rhs (and optionally lhs) terms.
 *
 *  \param flux_r    \ref Flux_Ref.
 *  \param volume    \ref Volume.
 *  \param s_store_i \ref Solver_Storage_Implicit.
 *  \param sim       \ref Simulation.
 */
typedef void (*compute_rlhs_fptr)
	(const struct Flux_Ref* flux_r,
	 struct Volume* volume,
	 struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim
	);

/// \brief Container for solver-related parameters.
struct S_Params {
	constructor_sol_vc_fptr constructor_sol_vc; ///< Pointer to the appropriate function.
	destructor_sol_vc_fptr  destructor_sol_vc;  ///< Pointer to the appropriate function.

	compute_rlhs_fptr compute_rlhs; ///< Pointer to the appropriate function.
};

/// \brief Container for the reference flux related parameters.
struct Flux_Ref {
	const struct const_Multiarray_d* fr;     ///< The reference flux.
	const struct const_Multiarray_d* dfr_ds; ///< The reference flux Jacobians with respect to the solution.
	const struct const_Multiarray_d* dfr_dg; ///< The reference flux Jacobians with respect to the solution gradients.
};

/** \brief Set the parameters of \ref S_Params.
 *  \return A statically allocated \ref S_Params container. */
static struct S_Params set_s_params
	(const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Flux_Ref container.
 *  \return See brief.
 *
 *  This function constructs the reference fluxes by multiplying the fluxes with the appropriate metric terms as
 *  specified by Zwanenburg et al. (eq. (B.3), \cite Zwanenburg2016). The memory layout of the reference flux terms
 *  (nodes,dim,eq) was chosen such that terms are grouped by dimension, allowing for differentiation operators to be
 *  applied efficiently; note that **this is not the same ordering as that used for the physical flux**. Please consult
 *  \ref compute_geometry_volume for the ordering of the metric terms if desired.
 */
static struct Flux_Ref* constructor_Flux_Ref
	(const struct const_Multiarray_d* m, ///< The metric terms.
	 const struct Flux* flux             ///< The physical \ref Flux.
	);

/// \brief Destructor for a \ref Flux_Ref container.
static void destructor_Flux_Ref
	(struct Flux_Ref* flux_ref ///< Standard.
	);

// Interface functions ********************************************************************************************** //

void compute_volume_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* s_store_i)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->elements->name == IL_ELEMENT_SOLVER_DG);

	struct S_Params s_params = set_s_params(sim);
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Volume*        vol   = (struct Volume*) curr;
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
//printf("%d\n",vol->index);

		// Compute the solution, gradients and xyz coordinates at the volume cubature nodes.
		flux_i->s   = s_params.constructor_sol_vc(vol,sim);
/// \todo Add functions computing gradients, xyz.
		flux_i->g   = NULL;
		flux_i->xyz = NULL;;
//print_Multiarray_d(s_vol->sol_coef);
//print_const_Multiarray_d(flux_i->s);

		// Compute the fluxes (and optionally their Jacobians) at the volume cubature nodes.
		struct Flux* flux = constructor_Flux(flux_i);
		s_params.destructor_sol_vc(flux_i->s);
//print_const_Multiarray_d(flux->f);
//print_const_Multiarray_d(flux->df_ds);

		// Compute the reference fluxes (and optionally their Jacobians) at the volume cubature nodes.
		struct Flux_Ref* flux_r = constructor_Flux_Ref(s_vol->metrics_vc,flux);
		destructor_Flux(flux);
//print_const_Multiarray_d(flux_r->dfr_ds);

		// Compute the rhs (and optionally the lhs) terms.
		s_params.compute_rlhs(flux_r,vol,s_store_i,sim);
		destructor_Flux_Ref(flux_r);
//EXIT_UNSUPPORTED;
	}
	destructor_Flux_Input(flux_i);
//EXIT_UNSUPPORTED;
}

struct Multiarray_Operator get_operator__tw1_vs_vc__rlhs_dg (const struct Volume* volume)
{
	struct Solver_Volume* s_volume = (struct Solver_Volume*) volume;
	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) volume->element;

	const int p      = s_volume->p_ref,
	          curved = volume->curved;
	const struct Multiarray_Operator tw1_vs_vc = {};
	set_MO_from_MO(&tw1_vs_vc,e->tw1_vs_vc[curved],1,(ptrdiff_t[]){0,0,p,p});

	return tw1_vs_vc;
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

/** \brief Constructor for a \ref const_Multiarray_d\* of reference flux from physical flux.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_flux_ref
	(const struct const_Multiarray_d* m, ///< Defined for \ref constructor_Flux_Ref.
	 const struct const_Multiarray_d* f  ///< The physical flux data.
	);

/// \brief Compute only the rhs term.
static void compute_rhs
	(const struct Flux_Ref* flux_r,             ///< Defined for \ref compute_rlhs_fptr.
	 struct Volume* volume,                     ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr.
	);

/// \brief Compute the rhs and lhs terms for 1st order equations only.
static void compute_rlhs_1
	(const struct Flux_Ref* flux_r,             ///< Defined for \ref compute_rlhs_fptr.
	 struct Volume* vol,                        ///< Defined for \ref compute_rlhs_fptr.
	 struct Solver_Storage_Implicit* s_store_i, ///< Defined for \ref compute_rlhs_fptr.
	 const struct Simulation* sim               ///< Defined for \ref compute_rlhs_fptr.
	);

static struct S_Params set_s_params (const struct Simulation* sim)
{
	struct S_Params s_params;

	if (!sim->collocated) {
		s_params.constructor_sol_vc = constructor_sol_vc_interp;
		s_params.destructor_sol_vc  = destructor_sol_vc_interp;
	} else {
		s_params.constructor_sol_vc = constructor_sol_vc_col;
		s_params.destructor_sol_vc  = destructor_sol_vc_col;
	}

	struct Test_Case* test_case = sim->test_case;
	switch (test_case->solver_method_curr) {
	case 'e':
		s_params.compute_rlhs = compute_rhs;
		break;
	case 'i':
		if (test_case->has_1st_order && !test_case->has_2nd_order)
			s_params.compute_rlhs = compute_rlhs_1;
		else if (!test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_2;
		else if (test_case->has_1st_order && test_case->has_2nd_order)
			EXIT_ADD_SUPPORT; //s_params.compute_rlhs = compute_rlhs_12;
		else
			EXIT_ERROR("Unsupported: %d %d\n",test_case->has_1st_order,test_case->has_2nd_order);
		break;
	default:
		EXIT_ERROR("Unsupported: %c\n",test_case->solver_method_curr);
		break;
	}

	return s_params;
}

static struct Flux_Ref* constructor_Flux_Ref (const struct const_Multiarray_d* m, const struct Flux* flux)
{
	assert(flux->f != NULL);
	assert(m->extents[0] == flux->f->extents[0]);

	struct Flux_Ref* flux_r = calloc(1,sizeof *flux_r); // returned

	flux_r->fr     = ( flux->f     ? constructor_flux_ref(m,flux->f)     : NULL );
	flux_r->dfr_ds = ( flux->df_ds ? constructor_flux_ref(m,flux->df_ds) : NULL );
	flux_r->dfr_dg = ( flux->df_dg ? constructor_flux_ref(m,flux->df_dg) : NULL );

	return flux_r;
}

static void destructor_Flux_Ref (struct Flux_Ref* flux_ref)
{
	if (flux_ref->fr)
		destructor_const_Multiarray_d(flux_ref->fr);
	if (flux_ref->dfr_ds)
		destructor_const_Multiarray_d(flux_ref->dfr_ds);
	if (flux_ref->dfr_dg)
		destructor_const_Multiarray_d(flux_ref->dfr_dg);
	free(flux_ref);
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_d* constructor_sol_vc_interp (struct Volume* volume, const struct Simulation* sim)
{
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	struct Solver_Volume* s_volume = (struct Solver_Volume*) volume;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) volume->element;

	const int p = s_volume->p_ref,
	          curved = volume->curved;
	const struct Operator* cv0_vs_vc = get_Multiarray_Operator(e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

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

static const struct const_Multiarray_d* constructor_flux_ref
	(const struct const_Multiarray_d* m, const struct const_Multiarray_d* f)
{
	assert(f->layout == 'C');

	const int order = f->order;
	ptrdiff_t extents[order];
	for (int i = 0; i < order; ++i) {
		if (i == 0)
			extents[i] = f->extents[i];
		else if (i == 1)
			extents[order-1] = f->extents[i];
		else
			extents[i-1] = f->extents[i];
	}

	struct Multiarray_d* fr = constructor_zero_Multiarray_d('C',order,extents); // returned

	const int n_n   = extents[0];
	const int d     = extents[order-1];
	const int n_col = compute_size(order,extents)/(n_n*d);

	int ind_f = 0;
	for (int col = 0; col < n_col; ++col) {
		for (int dim0 = 0; dim0 < d; ++dim0) {
			const int ind_fr = (ind_f+dim0*n_col);
			for (int dim1 = 0; dim1 < d; ++dim1) {
				const int ind_m  = dim0*d+dim1;
				const int ind_fp = (ind_f*d)+dim1;
				z_yxpz(n_n,get_col_const_Multiarray_d(ind_fp,f),
				           get_col_const_Multiarray_d(ind_m,m),
				           get_col_Multiarray_d(ind_fr,fr));
			}
		}
		++ind_f;
	}

	return (const struct const_Multiarray_d*) fr;
}

static void compute_rhs
	(const struct Flux_Ref* flux_r, struct Volume* volume, struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim)
{
	UNUSED(s_store_i);

	const struct Multiarray_Operator tw1_vs_vc = get_operator__tw1_vs_vc__rlhs_dg(volume);

	struct DG_Solver_Volume* dg_s_volume = (struct DG_Solver_Volume*) volume;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	const ptrdiff_t d = sim->d;
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vs_vc.data[dim],flux_r->fr,dg_s_volume->rhs,op_format,2,&dim,NULL);
}

static void compute_rlhs_1
	(const struct Flux_Ref* flux_r, struct Volume* vol, struct Solver_Storage_Implicit* s_store_i,
	 const struct Simulation* sim)
{
/// \todo Add special case for collocated.
// If collocation is enabled, note that the diagonal weight scaling must be added back in to recover the symmetry of the
// residual Jacobian. Add it just before adding the contribution to the petsc mat. Also add for face terms and RHS
// terms (volume, face, source or simply the complete rhs).
assert(sim->collocated == false); // Add support in future.
	const ptrdiff_t d = sim->d;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Solver_Volume* s_vol       = (struct Solver_Volume*) vol;
	struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) vol;

	const struct DG_Solver_Element* e = (const struct DG_Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	const struct Multiarray_Operator tw1_vs_vc;
	set_MO_from_MO(&tw1_vs_vc,e->tw1_vs_vc[curved],1,(ptrdiff_t[]){0,0,p,p});

	// rhs
	for (ptrdiff_t dim = 0; dim < d; ++dim)
		mm_NNC_Operator_Multiarray_d(1.0,1.0,tw1_vs_vc.data[dim],flux_r->fr,dg_s_vol->rhs,op_format,2,&dim,NULL);
//print_Multiarray_d(dg_s_vol->rhs);

	// lhs
	const struct Operator* cv0_vs_vc = get_Multiarray_Operator(e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const ptrdiff_t ext_0 = tw1_vs_vc.data[0]->op_std->ext_0,
	                ext_1 = tw1_vs_vc.data[0]->op_std->ext_1;

	struct Matrix_d* tw1_r = constructor_empty_Matrix_d('R',ext_0,ext_1);                    // destructed
	struct Matrix_d* lhs   = constructor_empty_Matrix_d('R',ext_0,cv0_vs_vc->op_std->ext_1); // destructed

	const struct const_Multiarray_d* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_d dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	const int n_eq = sim->test_case->n_eq,
	          n_vr = sim->test_case->n_var;
	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_d(tw1_r,0.0);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (double*)&dfr_ds_Ma->data[ind];
			mm_diag_d('R',1.0,1.0,tw1_vs_vc.data[dim]->op_std,(struct const_Vector_d*)&dfr_ds,tw1_r,false);
		}

		mm_d('N','N',1.0,0.0,(struct const_Matrix_d*)tw1_r,cv0_vs_vc->op_std,lhs);
//print_Matrix_d(lhs);

		set_petsc_Mat_row_col(s_store_i,s_vol,eq,s_vol,vr);
		add_to_petsc_Mat(s_store_i,(struct const_Matrix_d*)lhs);
	}}
	destructor_Matrix_d(tw1_r);
	destructor_Matrix_d(lhs);
}
