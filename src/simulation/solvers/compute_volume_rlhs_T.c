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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"


#include "def_templates_compute_volume_rlhs.h"

#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_flux.h"
#include "def_templates_math_functions.h"
#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref constructor_sol_vc_fptr using interpolation.
 *  \return See brief.. */
static const struct const_Multiarray_T* constructor_sol_vc_interp_T
	(const struct Solver_Volume_T* s_vol, ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/** \brief Version of \ref constructor_sol_vc_fptr using collocation.
 *  \return Standard. */
static const struct const_Multiarray_T* constructor_sol_vc_col_T
	(const struct Solver_Volume_T* s_vol, ///< See brief.
	 const struct Simulation* sim         ///< See brief.
	);

/// \brief Destructor for the return value of \ref constructor_sol_vc_interp_T.
static void destructor_sol_vc_interp_T
	(const struct const_Multiarray_T* sol_vc ///< To be destructed.
	);

/// \brief Destructor for the return value of \ref constructor_sol_vc_col_T.
static void destructor_sol_vc_col_T
	(const struct const_Multiarray_T* sol_vc ///< To be destructed.
	);

/** \brief Constructor for a \ref Flux_Ref_T container.
 *  \return See brief.
 *
 *  This function constructs the reference fluxes by multiplying the fluxes with the appropriate metric terms as
 *  specified by Zwanenburg et al. (eq. (B.3), \cite Zwanenburg2016).
 *
 *  The memory layout of the reference flux and reference flux Jacobian terms is that of the corresponding physical flux
 *  terms with the second extent (dim) moved to the last extent. For example:
 *  - df_ds: (nodes,dim,eq,var) -> dfr_ds: (nodes,eq,var,dim).
 *  This memory layout was chosen such that terms are grouped by dimension, allowing for differentiation operators to be
 *  applied efficiently; note that **this is not the same ordering as that used for the physical flux**. Please consult
 *  \ref compute_geometry_volume for the ordering of the metric terms if desired.
 */
static struct Flux_Ref_T* constructor_Flux_Ref_T
	(const struct const_Multiarray_R* m, ///< The metric terms.
	 const struct Flux_T* flux           ///< The physical \ref Flux_T.
	);

// Interface functions ********************************************************************************************** //

void set_S_Params_Volume_Structor_T (struct S_Params_Volume_Structor_T* spvs, const struct Simulation* sim)
{
	if (!sim->collocated) {
		spvs->constructor_sol_vc = constructor_sol_vc_interp_T;
		spvs->destructor_sol_vc  = destructor_sol_vc_interp_T;
	} else {
		spvs->constructor_sol_vc = constructor_sol_vc_col_T;
		spvs->destructor_sol_vc  = destructor_sol_vc_col_T;
	}
}

struct Flux_Ref_T* constructor_Flux_Ref_vol_T
	(const struct S_Params_Volume_Structor_T* spvs, struct Flux_Input_T* flux_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	// Compute the solution, gradients and xyz coordinates at the volume cubature nodes.
	flux_i->s   = spvs->constructor_sol_vc(s_vol,sim);
/// \todo Add functions computing gradients, xyz.
	flux_i->g   = NULL;
	flux_i->xyz = NULL;
//print_Multiarray_T(s_vol->sol_coef);
//print_const_Multiarray_T(flux_i->s);

	// Compute the fluxes (and optionally their Jacobians) at the volume cubature nodes.
	struct Flux_T* flux = constructor_Flux_T(flux_i); // destructed
	spvs->destructor_sol_vc(flux_i->s);
//print_const_Multiarray_T(flux->f);
//print_const_Multiarray_T(flux->df_ds);

	// Compute the reference fluxes (and optionally their Jacobians) at the volume cubature nodes.
	struct Flux_Ref_T* flux_r = constructor_Flux_Ref_T(s_vol->metrics_vc,flux);
//print_const_Multiarray_T(flux_r->fr);
//print_const_Multiarray_T(flux_r->dfr_ds);
	destructor_Flux_T(flux);
//EXIT_UNSUPPORTED;

	return flux_r;
}

void destructor_Flux_Ref_T (struct Flux_Ref_T* flux_ref)
{
/// \todo Make a conditional destructor.
	if (flux_ref->fr)
		destructor_const_Multiarray_T(flux_ref->fr);
	if (flux_ref->dfr_ds)
		destructor_const_Multiarray_T(flux_ref->dfr_ds);
	if (flux_ref->dfr_dg)
		destructor_const_Multiarray_T(flux_ref->dfr_dg);
	if (flux_ref->d2fr_ds2)
		destructor_const_Multiarray_T(flux_ref->d2fr_ds2);
	free(flux_ref);
}

struct Matrix_T* constructor_lhs_v_1_T
	(const struct Flux_Ref_T* flux_r, const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	const struct Multiarray_Operator tw1_vt_vc = get_operator__tw1_vt_vc_T(s_vol);
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc_T(s_vol);

	const ptrdiff_t ext_0 = tw1_vt_vc.data[0]->op_std->ext_0,
	                ext_1 = tw1_vt_vc.data[0]->op_std->ext_1;
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_eq = test_case->n_eq,
	          n_vr = test_case->n_var;

	struct Matrix_T* tw1_r = constructor_empty_Matrix_T('R',ext_0,ext_1);                         // destructed
	struct Matrix_T* lhs_l = constructor_empty_Matrix_T('R',ext_0,cv0_vs_vc->op_std->ext_1);      // destructed
	struct Matrix_T* lhs   = constructor_empty_Matrix_T('R',n_eq*lhs_l->ext_0,n_vr*lhs_l->ext_1); // returned

	const struct const_Multiarray_T* dfr_ds_Ma = flux_r->dfr_ds;
	struct Vector_T dfr_ds = { .ext_0 = dfr_ds_Ma->extents[0], .owns_data = false, .data = NULL, };

	for (int vr = 0; vr < n_vr; ++vr) {
	for (int eq = 0; eq < n_eq; ++eq) {
		set_to_value_Matrix_T(tw1_r,0.0);
		for (int dim = 0; dim < DIM; ++dim) {
			const ptrdiff_t ind =
				compute_index_sub_container(dfr_ds_Ma->order,1,dfr_ds_Ma->extents,(ptrdiff_t[]){eq,vr,dim});
			dfr_ds.data = (Type*)&dfr_ds_Ma->data[ind];
			mm_diag_T('R',1.0,1.0,tw1_vt_vc.data[dim]->op_std,(struct const_Vector_T*)&dfr_ds,tw1_r,false);
		}

		mm_TRT('N','N',1.0,0.0,(struct const_Matrix_T*)tw1_r,cv0_vs_vc->op_std,lhs_l);
//printf("%d %d\n",vr,eq);
//print_Matrix_T(lhs_l);
		set_block_Matrix_T(lhs,(struct const_Matrix_T*)lhs_l,eq*lhs_l->ext_0,vr*lhs_l->ext_1,'i');
	}}
	destructor_Matrix_T(tw1_r);
	destructor_Matrix_T(lhs_l);

	return lhs;
}

const struct Operator* get_operator__cv0_vs_vc_T (const struct Solver_Volume_T* s_vol)
{
	const struct Volume* vol       = (struct Volume*) s_vol;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int p = s_vol->p_ref,
	          curved = vol->curved;
	return get_Multiarray_Operator(e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});
}

struct Multiarray_Operator get_operator__tw1_vt_vc_T (const struct Solver_Volume_T* s_vol)
{
	const struct Volume* vol       = (struct Volume*) s_vol;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int p      = s_vol->p_ref,
	          curved = vol->curved;
	struct Multiarray_Operator tw1_vt_vc = set_MO_from_MO(e->tw1_vt_vc[curved],1,(ptrdiff_t[]){0,0,p,p});

	return tw1_vt_vc;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for a \ref const_Multiarray_T\* of reference flux from physical flux.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_flux_ref_T
	(const struct const_Multiarray_R* m, ///< Defined for \ref constructor_Flux_Ref_T.
	 const struct const_Multiarray_T* f  ///< The physical flux data.
	);

static const struct const_Multiarray_T* constructor_sol_vc_interp_T
	(const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	const struct Operator* cv0_vs_vc = get_operator__cv0_vs_vc_T(s_vol);

	const struct const_Multiarray_T* s_coef = (const struct const_Multiarray_T*) s_vol->sol_coef;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';
	return constructor_mm_NN1_Operator_const_Multiarray_T(cv0_vs_vc,s_coef,'C',op_format,s_coef->order,NULL);
}

static const struct const_Multiarray_T* constructor_sol_vc_col_T
	(const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	UNUSED(sim);
	return (const struct const_Multiarray_T*) s_vol->sol_coef;
}

static void destructor_sol_vc_interp_T (const struct const_Multiarray_T* sol_vc)
{
	destructor_const_Multiarray_T(sol_vc);
}

static void destructor_sol_vc_col_T (const struct const_Multiarray_T* sol_vc)
{
	UNUSED(sol_vc);
}

static struct Flux_Ref_T* constructor_Flux_Ref_T (const struct const_Multiarray_R* m, const struct Flux_T* flux)
{
	assert(flux->f != NULL);
	assert(m->extents[0] == flux->f->extents[0]);

	struct Flux_Ref_T* flux_r = calloc(1,sizeof *flux_r); // returned

	flux_r->fr       = ( flux->f       ? constructor_flux_ref_T(m,flux->f)       : NULL );
	flux_r->dfr_ds   = ( flux->df_ds   ? constructor_flux_ref_T(m,flux->df_ds)   : NULL );
	flux_r->dfr_dg   = ( flux->df_dg   ? constructor_flux_ref_T(m,flux->df_dg)   : NULL );
	flux_r->d2fr_ds2 = ( flux->d2f_ds2 ? constructor_flux_ref_T(m,flux->d2f_ds2) : NULL );

	return flux_r;
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_T* constructor_flux_ref_T
	(const struct const_Multiarray_R* m, const struct const_Multiarray_T* f)
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

	struct Multiarray_T* fr = constructor_zero_Multiarray_T('C',order,extents); // returned

	const int n_n   = (int)extents[0];
	const int n_col = (int)compute_size(order,extents)/(n_n*DIM);
	assert(extents[order-1] == DIM);

	int ind_f = 0;
	for (int col = 0; col < n_col; ++col) {
		for (int dim0 = 0; dim0 < DIM; ++dim0) {
			const int ind_fr = (ind_f+dim0*n_col);
			for (int dim1 = 0; dim1 < DIM; ++dim1) {
				const int ind_m  = dim0*DIM+dim1,
				          ind_fp = (ind_f*DIM)+dim1;
				z_yxpz_RTT(n_n,get_col_const_Multiarray_R(ind_m,m),
				               get_col_const_Multiarray_T(ind_fp,f),
				               get_col_Multiarray_T(ind_fr,fr));
			}
		}
		++ind_f;
	}

	return (const struct const_Multiarray_T*) fr;
}
