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

#include "solution.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "matrix.h"

#include "computational_elements.h"
#include "face_solver.h"
#include "element.h"
#include "element_solution.h"
#include "volume.h"
#include "volume_solver.h"

#include "flux.h"
#include "geometry.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set up the initial \ref Solver_Volume::sol_coef and \ref Solver_Volume::grad_coef.
static void set_initial_v_sg_coef
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set up the initial \ref Solver_Face::nf_coef.
static void set_initial_f_nf_coef
	(struct Simulation* sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_d* constructor_const_sol_invalid
	(const struct const_Multiarray_d* xyz, const struct Simulation* sim)
{
	UNUSED(xyz);
	UNUSED(sim);
	EXIT_ERROR("Should not be entering here.\n");
	return NULL;
}

void set_initial_solution (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);

	constructor_derived_Elements(sim,IL_ELEMENT_SOLUTION);

	switch (sim->method) {
	case METHOD_DG:
		set_initial_v_sg_coef(sim);
		break;
	case METHOD_DPG:
		set_initial_v_sg_coef(sim);
		set_initial_f_nf_coef(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	destructor_derived_Elements(sim,IL_ELEMENT);
}

void set_sg_do_nothing (const struct Simulation* sim, struct Solution_Container sol_cont)
{
	UNUSED(sim);
	UNUSED(sol_cont);
	return;
}

const struct const_Multiarray_d* constructor_xyz_v
	(const struct Simulation* sim, struct Solver_Volume* volume, const char node_kind)
{
UNUSED(sim);
	assert((node_kind == 's') || (node_kind == 'c'));

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) volume;
	const struct Solution_Element* element = (struct Solution_Element*) base_volume->element;

	const int d = ((struct const_Element*)element)->d;

	const int curved = base_volume->curved,
	          p_o    = volume->p_ref,
	          p_i    = ( curved ? p_o : 1 );

	const struct Operator* cv0_vg_vX = NULL;
	if (node_kind == 's')
		cv0_vg_vX = get_Multiarray_Operator(element->cv0_vg_vs[curved],(ptrdiff_t[]){0,0,p_o,p_i});
	else if (node_kind == 'c')
		cv0_vg_vX = get_Multiarray_Operator(element->cv0_vg_vc[curved],(ptrdiff_t[]){0,0,p_o,p_i});

	const int n_vs = cv0_vg_vX->op_std->ext_0;

	const struct const_Multiarray_d*const geom_coef = volume->geom_coef;
	struct Multiarray_d* xyz_v = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_vs,d}); // returned
	mm_NN1C_Operator_Multiarray_d(cv0_vg_vX,geom_coef,xyz_v,op_format,geom_coef->order,NULL,NULL);

	return (const struct const_Multiarray_d*) xyz_v;
}

const struct const_Multiarray_d* constructor_xyz_f
	(const struct Simulation* sim, struct Solver_Face* s_face, const char node_kind)
{
	UNUSED(sim);
	assert((node_kind == 'f') || (node_kind == 't'));

	struct Face* face  = (struct Face*) s_face;
	struct Volume* vol = (struct Volume*) face->neigh_info[0].volume;

	const struct Solution_Element* e = (struct Solution_Element*) vol->element;

	const int ind_lf   = face->neigh_info[0].ind_lf,
	          ind_href = face->neigh_info[0].ind_href,
	          curved   = vol->curved;
	const int p_v = ( curved ? ((struct Solver_Volume*)vol)->p_ref : 1 ),
	          p_f = s_face->p_ref;

	const struct Operator* cv0_vg_fX = NULL;
	if (node_kind == 'f')
		cv0_vg_fX = get_Multiarray_Operator(e->cv0_vg_ff[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
	else if (node_kind == 't')
		EXIT_ADD_SUPPORT;

	const int d    = ((struct const_Element*)e)->d,
	          n_fs = cv0_vg_fX->op_std->ext_0;

	struct Multiarray_d* xyz_f = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_fs,d}); // returned

	const char op_format = 'd';

	const struct const_Multiarray_d*const geom_coef = ((struct Solver_Volume*)vol)->geom_coef;
	mm_NN1C_Operator_Multiarray_d(cv0_vg_fX,geom_coef,xyz_f,op_format,geom_coef->order,NULL,NULL);

	return (const struct const_Multiarray_d*) xyz_f;
}

void compute_coef_from_val_vs
	(const struct Solver_Volume* s_vol, const struct const_Multiarray_d* sol_val, struct Multiarray_d* sol_coef)
{
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) s_vol;
	const struct Solution_Element* element = (struct Solution_Element*) base_volume->element;

	const int p_ref = s_vol->p_ref;

	const struct Operator* vc0_vs_vs = get_Multiarray_Operator(element->vc0_vs_vs,(ptrdiff_t[]){0,0,p_ref,p_ref});

	resize_Multiarray_d(sol_coef,sol_val->order,sol_val->extents);
	mm_NN1C_Operator_Multiarray_d(vc0_vs_vs,sol_val,sol_coef,op_format,sol_coef->order,NULL,NULL);
}

struct Multiarray_d* constructor_sol_v
	(const struct Simulation* sim, struct Solver_Volume* volume, const char node_kind)
{
UNUSED(sim);
	assert((node_kind == 's') || (node_kind == 'c'));

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* base_volume = (struct Volume*) volume;
	const struct Solution_Element* element = (struct Solution_Element*) base_volume->element;

	const int curved = base_volume->curved,
	          p      = volume->p_ref;

	const struct Operator* cv0_vs_vX = NULL;
	if (node_kind == 's')
		cv0_vs_vX = get_Multiarray_Operator(element->cv0_vs_vs,(ptrdiff_t[]){0,0,p,p});
	else if (node_kind == 'c')
		cv0_vs_vX = get_Multiarray_Operator(element->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const struct const_Multiarray_d*const sol_coef = (struct const_Multiarray_d*) volume->sol_coef;

	const ptrdiff_t ext_0 = cv0_vs_vX->op_std->ext_0,
	                ext_1 = sol_coef->extents[1];

	struct Multiarray_d* sol_v = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,ext_1}); // returned
	mm_NN1C_Operator_Multiarray_d(cv0_vs_vX,sol_coef,sol_v,op_format,sol_coef->order,NULL,NULL);

	return sol_v;
}

void compute_source_do_nothing (const struct Simulation* sim, struct Solver_Volume* volume)
{
	UNUSED(sim);
	UNUSED(volume);
	return;
}

void update_Solution_Container_sol (struct Solution_Container*const sol_cont, struct Multiarray_d*const sol)
{
	const char cv_type = sol_cont->cv_type;

	if (cv_type == 'v') {
		assert(sol_cont->sol->data != NULL);

		sol_cont->sol->extents[0] = sol->extents[0];
		sol_cont->sol->extents[1] = sol->extents[1];
		free(sol_cont->sol->data);
		sol_cont->sol->data = sol->data;

		sol->owns_data = false;
	} else if (cv_type == 'c') {
		assert(sol_cont->node_kind == 's');
		compute_coef_from_val_vs(sol_cont->volume,(struct const_Multiarray_d*)sol,sol_cont->sol);
	} else {
		EXIT_ERROR("Unsupported: %c\n",cv_type);
	}
}


// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the multiarray of normal terms at the face flux nodes.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_normals_ff
	(const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

/** \brief Constructor for the multiarray of normal flux terms.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_nf
	(const struct const_Multiarray_d* flux,   ///< The flux data.
	 const struct const_Multiarray_d* normals ///< The unit normals data.
	);

/// \brief Compute the coefficients associated with the values of the face normal flux.
void compute_coef_from_val_ff
	(const struct Solver_Face* s_face,       ///< \ref Solver_Face.
	 const struct const_Multiarray_d* f_val, ///< The flux values.
	 struct Multiarray_d* f_coef             ///< To hold the flux coefficients.
	);

static void set_initial_v_sg_coef (struct Simulation* sim)
{
	struct Solution_Container sol_cont =
		{ .ce_type = 'v', .cv_type = 'c', .node_kind = 's', .volume = NULL, .face = NULL, .sol = NULL, };
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* volume = (struct Solver_Volume*) curr;

		sol_cont.volume = volume;
		sol_cont.sol = volume->sol_coef;
		sim->test_case->set_sol(sim,sol_cont);

		sol_cont.sol = volume->grad_coef;
		sim->test_case->set_grad(sim,sol_cont);
	}
}

static void set_initial_f_nf_coef (struct Simulation* sim)
{
	struct Multiarray_d* sol_fs = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}); // destructed

	sim->test_case->solver_method_curr = 'e';
	struct Flux_Input* flux_i = constructor_Flux_Input(sim); // destructed
	sim->test_case->solver_method_curr = '0';

	struct Solution_Container sol_cont =
		{ .ce_type = 'f', .cv_type = 'v', .node_kind = 'f', .volume = NULL, .face = NULL, .sol = NULL, };
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Solver_Face* s_face = (struct Solver_Face*) curr;

		sol_cont.face = s_face;
		sol_cont.sol  = sol_fs;
		sim->test_case->set_sol(sim,sol_cont);

		flux_i->s = (struct const_Multiarray_d*)sol_fs;

		struct Flux* flux = constructor_Flux(flux_i);

		const struct const_Multiarray_d* normals_ff = constructor_normals_ff(s_face); // destructed

		const struct const_Multiarray_d* nff = constructor_nf(flux->f,normals_ff); // destructed
		destructor_Flux(flux);
		destructor_const_Multiarray_d(normals_ff);

		compute_coef_from_val_ff(s_face,nff,s_face->nf_coef);
		destructor_const_Multiarray_d(nff);
	}
	destructor_Multiarray_d(sol_fs);

	destructor_Flux_Input(flux_i);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the multiarray of metric terms at the face flux nodes.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_metrics_ff
	(const struct Solver_Face* s_face ///< \ref Solver_Face.
	);

static const struct const_Multiarray_d* constructor_normals_ff (const struct Solver_Face* s_face)
{
	const struct const_Multiarray_d* metrics_ff = constructor_metrics_ff(s_face); // destructed

	const struct Face* face       = (struct Face*) s_face;
	const struct Volume* vol      = face->neigh_info[0].volume;
	const struct const_Element* e = (struct const_Element*) vol->element;

	struct Multiarray_d* normals_ff = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}); // returned
	compute_unit_normals(face->neigh_info[0].ind_lf,e->normals,metrics_ff,normals_ff);
	destructor_const_Multiarray_d(metrics_ff);

	transpose_Multiarray_d(normals_ff,true); // Convert to column major.

	return (struct const_Multiarray_d*) normals_ff;
}

static const struct const_Multiarray_d* constructor_nf
	(const struct const_Multiarray_d* flux, const struct const_Multiarray_d* normals)
{
	assert(flux->layout == 'C');
	assert(flux->extents[0] == normals->extents[0]);
	assert(flux->extents[1] == normals->extents[1]);
	assert(flux->layout == normals->layout);

	const int n_n  = flux->extents[0],
	          d    = flux->extents[1],
	          n_eq = flux->extents[2];

	struct Multiarray_d* nf = constructor_zero_Multiarray_d('C',2,(ptrdiff_t[]){n_n,n_eq}); // returned

	for (int eq = 0; eq < n_eq; ++eq) {
		double*const data_nf = get_col_Multiarray_d(eq,nf);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind_f = compute_index_sub_container(flux->order,1,flux->extents,(ptrdiff_t[]){dim,eq});
			const double*const data_f = get_col_const_Multiarray_d(ind_f,flux),
			            *const data_n = get_col_const_Multiarray_d(dim,normals);
			for (int n = 0; n < n_n; ++n)
				data_nf[n] += data_f[n]*data_n[n];
		}
	}

	return (struct const_Multiarray_d*) nf;
}

void compute_coef_from_val_ff
	(const struct Solver_Face* s_face, const struct const_Multiarray_d* f_val, struct Multiarray_d* f_coef)
{
	const struct Face* face          = (struct Face*) s_face;
	const struct Solution_Element* e = (struct Solution_Element*) face->neigh_info[0].volume->element;

	const int ind_e = get_face_element_index(face),
	          p     = s_face->p_ref;
	const struct Operator* vc0_ff_ff = get_Multiarray_Operator(e->vc0_ff_ff,(ptrdiff_t[]){ind_e,ind_e,0,0,p,p});

	resize_Multiarray_d(f_coef,f_val->order,f_val->extents);
	mm_NN1C_Operator_Multiarray_d(vc0_ff_ff,f_val,f_coef,'d',f_coef->order,NULL,NULL);
}

// Level 2 ********************************************************************************************************** //

static const struct const_Multiarray_d* constructor_metrics_ff (const struct Solver_Face* s_face)
{
	const struct Face* face           = (struct Face*) s_face;
	const struct Volume* vol          = face->neigh_info[0].volume;
	const struct Solver_Volume* s_vol = (struct Solver_Volume*) vol;

	const struct Solution_Element* e = (struct Solution_Element*) vol->element;

	const int ind_lf = face->neigh_info[0].ind_lf,
	          cv     = vol->curved,
	          pv     = ( cv ? s_vol->p_ref : 1),
	          pf     = s_face->p_ref;
	const struct Operator* vv0_vm_ff = get_Multiarray_Operator(e->vv0_vm_ff[cv],(ptrdiff_t[]){ind_lf,0,0,pf,pv});

	const struct const_Multiarray_d* m_vm = s_vol->metrics_vm;
	return constructor_mm_NN1_Operator_const_Multiarray_d(vv0_vm_ff,m_vm,'C','d',m_vm->order,NULL);
}
