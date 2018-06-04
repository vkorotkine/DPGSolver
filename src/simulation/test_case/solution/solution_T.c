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


#include "def_templates_solution.h"

#include "def_templates_multiarray.h"

#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

#include "def_templates_flux.h"
#include "def_templates_operators.h"
#include "def_templates_solve.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Set up the initial \ref Solver_Volume_T::sol_coef and \ref Solver_Volume_T::grad_coef.
static void set_initial_v_sg_coef
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Set up the initial \ref Solver_Face_T::nf_coef.
 *  \note In order to allow for the general imposition of boundary conditions which are dependent on the internal
 *        solution, it was decided to set normal flux unknowns on boundary faces based on the solution. This means that
 *        there are no `nf_coef` terms on boundary faces. */
static void set_initial_f_nf_coef
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set up the exact \ref Solver_Face_T::nf_fc if possible.
static void set_exact_f_nf_fc
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Set up the initial \ref Solver_Volume_T::test_s_coef and \todo ref Solver_Volume_T::test_g_coef.
static void set_initial_v_test_sg_coef
	(struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Get the pointer to the appropriate \ref Solver_Element::cv0_vg_vc operator.
 *  \return See brief. */
static const struct Operator* get_operator__cv0_vg_vc
	(const struct Solver_Volume_T* s_vol ///< The current volume.
	);

/** \brief Contructor for a \ref const_Multiarray_T\* holding the xyz coordinates at volume nodes of input kind.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_v
	(const struct Simulation* sim,  ///< \ref Simulation.
	 struct Solver_Volume_T* s_vol, ///< \ref Solver_Volume_T.
	 const char node_kind,          ///< The kind of node. Options: 's'olution, 'c'ubature.
	 const bool using_restart       ///< Flag for whether a restart solution is being used.
	);

/** \brief Contructor for a \ref const_Multiarray_T\* holding the xyz coordinates at face nodes of input kind.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_f
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Face_T* s_face, ///< \ref Solver_Face_T.
	 const char node_kind,         ///< The kind of node. Options: 'f'lux, 't'race.
	 const bool using_restart      ///< Flag for whether a restart solution is being used.
	);

/// \brief Compute the coefficients associated with the values of the volume solution.
static void compute_coef_from_val_vs
	(const struct Solver_Volume_T* s_vol,        ///< \ref Solver_Volume_T.
	 const struct const_Multiarray_T* sol_val, ///< The solution values.
	 struct Multiarray_T* sol_coef             ///< To hold the solution coefficients.
	);

/// \brief Compute the coefficients associated with the values of the volume gradient.
static void compute_coef_from_val_vg
	(const struct Solver_Volume_T* s_vol,       ///< \ref Solver_Volume_T.
	 const struct const_Multiarray_T* grad_val, ///< The gradient values.
	 struct Multiarray_T* grad_coef             ///< To hold the gradient coefficients.
	);

/** \brief Constructor for the multiarray of normal flux data evaluated at the desired face nodes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_nf
	(struct Solver_Face_T* s_face,     ///< \ref Solver_Face_T.
	 const char node_kind,             ///< The type of face node at which to obtain the normal flux values.
	 struct Flux_Input_T*const flux_i, ///< \ref Flux_Input_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Compute the coefficients associated with the values of the face normal flux.
static void compute_coef_from_val_ff
	(const struct Solver_Face_T* s_face,     ///< \ref Solver_Face_T.
	 const struct const_Multiarray_T* f_val, ///< The flux values.
	 struct Multiarray_T* f_coef             ///< To hold the flux coefficients.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_T* constructor_const_sol_invalid_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	UNUSED(xyz);
	UNUSED(sim);
	EXIT_ERROR("Should not be entering here.\n");
	return NULL;
}

void set_initial_solution_T (struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	switch (sim->method) {
	case METHOD_DG:
		set_initial_v_sg_coef(sim);
		break;
	case METHOD_OPG:
		set_initial_v_test_sg_coef(sim);
		// fallthrough
	case METHOD_DPG:
		set_initial_v_sg_coef(sim);
		set_initial_f_nf_coef(sim);
		set_exact_f_nf_fc(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
}

void set_sg_do_nothing_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(list_is_derived_from("solver",'e',sim));
	UNUSED(sol_cont);
	return;
}

void set_sg_zero_T (const struct Simulation*const sim, struct Solution_Container_T sol_cont)
{
	EXIT_ADD_SUPPORT; UNUSED(sim); UNUSED(sol_cont);
}

const struct const_Multiarray_R* constructor_xyz_sol_T
	(const struct Simulation* sim, const struct Solution_Container_T* sol_cont)
{
	assert(list_is_derived_from("solver",'e',sim));

	const char ce_type   = sol_cont->ce_type,
	           node_kind = sol_cont->node_kind;
	const bool using_restart = sol_cont->using_restart;

	const struct const_Multiarray_R* xyz = NULL;
	if (ce_type == 'v')
		xyz = constructor_xyz_v(sim,sol_cont->volume,node_kind,using_restart); // returned
	else if (ce_type == 'f')
		xyz = constructor_xyz_f(sim,sol_cont->face,node_kind,using_restart); // returned
	else
		EXIT_ERROR("Unsupported: %c\n",ce_type);

	return xyz;
}

struct Multiarray_T* constructor_sol_v_T
	(const struct Simulation* sim, struct Solver_Volume_T* s_vol, const char node_kind)
{
	assert(list_is_derived_from("solver",'e',sim));
	assert((node_kind == 's') || (node_kind == 'c'));

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int curved = vol->curved,
	          p      = s_vol->p_ref;

	const struct Operator* cv0_vs_vX = NULL;
	if (node_kind == 's')
		cv0_vs_vX = get_Multiarray_Operator(s_e->cv0_vs_vs,(ptrdiff_t[]){0,0,p,p});
	else if (node_kind == 'c')
		cv0_vs_vX = get_Multiarray_Operator(s_e->cv0_vs_vc[curved],(ptrdiff_t[]){0,0,p,p});

	const struct const_Multiarray_T*const sol_coef = (struct const_Multiarray_T*) s_vol->sol_coef;

	const ptrdiff_t ext_0 = cv0_vs_vX->op_std->ext_0,
	                ext_1 = sol_coef->extents[1];

	struct Multiarray_T* sol_v = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){ext_0,ext_1}); // returned
	mm_NN1C_Operator_Multiarray_T(cv0_vs_vX,sol_coef,sol_v,op_format,sol_coef->order,NULL,NULL);

	return sol_v;
}

void compute_source_rhs_do_nothing_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	UNUSED(sim);
	UNUSED(s_vol);
	UNUSED(rhs);
	return;
}

void add_to_flux_imbalance_source_do_nothing_T
	(const struct Simulation* sim, const struct Solver_Volume_T* s_vol, struct Multiarray_T* rhs)
{
	UNUSED(sim);
	UNUSED(s_vol);
	UNUSED(rhs);
	return;
}

void update_Solution_Container_sol_T
	(struct Solution_Container_T*const sol_cont, struct Multiarray_T*const sol, const struct Simulation*const sim)
{
	const int order_exp = 2;
	assert(list_is_derived_from("solver",'e',sim));
	assert(sol_cont->sol->order == sol->order);
	assert(sol_cont->sol->order == order_exp);

	const char cv_type = sol_cont->cv_type;

	if (cv_type == 'v') {
		assert(sol_cont->sol->data != NULL);

		for (int i = 0; i < order_exp; ++i)
			sol_cont->sol->extents[i] = sol->extents[i];
		free(sol_cont->sol->data);
		sol_cont->sol->data = sol->data;

		sol->owns_data = false;
	} else if (cv_type == 'c') {
		assert(sol_cont->node_kind == 's');
		compute_coef_from_val_vs(sol_cont->volume,(struct const_Multiarray_T*)sol,sol_cont->sol);
	} else {
		EXIT_ERROR("Unsupported: %c\n",cv_type);
	}
}

void update_Solution_Container_grad_T
	(struct Solution_Container_T*const sol_cont, struct Multiarray_T*const grad, const struct Simulation*const sim)
{
	const int order_exp = 3;
	assert(list_is_derived_from("solver",'e',sim));
	assert(sol_cont->sol->order == grad->order);
	assert(sol_cont->sol->order == order_exp);

	const char cv_type = sol_cont->cv_type;

	if (cv_type == 'v') {
		assert(sol_cont->sol->data != NULL);

		for (int i = 0; i < order_exp; ++i)
			sol_cont->sol->extents[i] = grad->extents[i];
		free(sol_cont->sol->data);
		sol_cont->sol->data = grad->data;

		grad->owns_data = false;
	} else if (cv_type == 'c') {
		assert(sol_cont->node_kind == 'r');
		compute_coef_from_val_vg(sol_cont->volume,(struct const_Multiarray_T*)grad,sol_cont->sol);
	} else {
		EXIT_ERROR("Unsupported: %c\n",cv_type);
	}
}

const struct const_Multiarray_R* constructor_xyz_vc_interp_T
	(const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	assert(list_is_derived_from("solver",'e',sim));
	const struct Operator* cv0_vg_vc = get_operator__cv0_vg_vc(s_vol);

	const struct const_Multiarray_R* geom_coef = s_vol->geom_coef;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';
	return constructor_mm_NN1_Operator_const_Multiarray_R(cv0_vg_vc,geom_coef,'C',op_format,geom_coef->order,NULL);
}

void constructor_Solver_Face__nf_coef_T
	(struct Solver_Face_T*const s_face, struct Flux_Input_T*const flux_i, const struct Simulation*const sim)
{
	const struct const_Multiarray_T*const nf = constructor_nf(s_face,'f',flux_i,sim); // destructed
	compute_coef_from_val_ff(s_face,nf,s_face->nf_coef);
	destructor_const_Multiarray_T(nf);
}

void set_grad_zero_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(sol_cont.node_kind != 's');
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* grad = (struct Multiarray_T*) constructor_const_grad_zero_T(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_grad_T(&sol_cont,grad,sim);
	destructor_Multiarray_T(grad);
}

const struct const_Multiarray_T* constructor_const_grad_zero_T
	(const struct const_Multiarray_R*const xyz, const struct Simulation*const sim)
{
	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	return (struct const_Multiarray_T*) constructor_zero_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_var,DIM});
}

void set_to_zero_residual_T (const struct Simulation*const sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		struct Multiarray_T*const s_coef = s_vol->sol_coef;

		destructor_conditional_Multiarray_T(s_vol->rhs_0);
		s_vol->rhs_0 = constructor_zero_Multiarray_T(s_coef->layout,s_coef->order,s_coef->extents); // keep
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the multiarray of normal terms at the face flux nodes.
 *  \return See brief. */
/// \todo move down one level.
static const struct const_Multiarray_R* constructor_normals_ff
	(const struct Solver_Face_T* s_face ///< \ref Solver_Face_T.
	);

/** \brief Constructor for the multiarray of normal flux terms as the dot product of the input normals and fluxes.
 *  \return See brief. */
static const struct const_Multiarray_T* constructor_n_dot_f
	(const struct const_Multiarray_T* flux,   ///< The flux data.
	 const struct const_Multiarray_R* normals ///< The unit normals data.
	);

static void set_initial_v_sg_coef (struct Simulation* sim)
{
	struct Solution_Container_T sol_cont =
		{ .ce_type = 'v', .cv_type = 'c', .node_kind = 0, .volume = NULL, .face = NULL, .sol = NULL, };
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;

		struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
		sol_cont.volume = s_vol;

		sol_cont.node_kind = 's';
		sol_cont.sol       = s_vol->sol_coef;
		test_case->set_sol_start(sim,sol_cont);

		sol_cont.node_kind = 'r';
		sol_cont.sol       = s_vol->grad_coef;
		test_case->set_grad(sim,sol_cont);
	}
}

static void set_initial_f_nf_coef (struct Simulation* sim)
{
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T_e(sim); // destructed
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		/// Boundary faces do not have normal flux unknowns.
		if (face->boundary)
			continue;

		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;
		constructor_Solver_Face__nf_coef_T(s_face,flux_i,sim);
	}
	destructor_Flux_Input_T(flux_i);
}

static void set_exact_f_nf_fc (struct Simulation* sim)
{
	struct Flux_Input_T* flux_i = constructor_Flux_Input_T_e(sim); // destructed
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		struct Face* face = (struct Face*) curr;
		/// Exact normal flux currently only being used on boundary faces.
		if (!face->boundary)
			continue;

		struct Solver_Face_T* s_face = (struct Solver_Face_T*) curr;
		s_face->nf_fc = constructor_nf(s_face,'c',flux_i,sim); // keep
	}
	destructor_Flux_Input_T(flux_i);
}

static void set_initial_v_test_sg_coef (struct Simulation*const sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) curr;
		const struct Operator*const tw0_vt_vc = get_operator__tw0_vt_vc_T(s_vol);

		const ptrdiff_t ext_0 = tw0_vt_vc->op_std->ext_0,
		                extents[] = { ext_0, n_var, };
		resize_Multiarray_T(s_vol->test_s_coef,s_vol->test_s_coef->order,extents);
		set_to_value_Multiarray_T(s_vol->test_s_coef,0.0);
	}
}

static const struct Operator* get_operator__cv0_vg_vc (const struct Solver_Volume_T* s_vol)
{
	const struct Volume* vol       = (struct Volume*) s_vol;
	const struct Solver_Element* e = (struct Solver_Element*) vol->element;

	const int p_o = s_vol->p_ref,
	          curved = vol->curved,
	          p_i = ( curved ? p_o : 1 );
	return get_Multiarray_Operator(e->cv0_vg_vc[curved],(ptrdiff_t[]){0,0,p_o,p_i});
}

static const struct const_Multiarray_R* constructor_xyz_v
	(const struct Simulation* sim, struct Solver_Volume_T* s_vol, const char node_kind, const bool using_restart)
{
UNUSED(sim);
	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int d = ((struct const_Element*)s_e)->d;

	const int curved = ( !using_restart ? vol->curved : false ),
	          p_o    = s_vol->p_ref,
	          p_i    = ( !curved ? 1 : p_o );

	const ptrdiff_t op_inds[] = {0,0,p_o,p_i};
	const struct Operator* cv0_vg_vX = NULL;
	switch (node_kind) {
		case 's': cv0_vg_vX = get_Multiarray_Operator(s_e->cv0_vg_vs[curved],op_inds); break;
		case 'r': cv0_vg_vX = get_Multiarray_Operator(s_e->cv0_vg_vr[curved],op_inds); break;
		case 'c': cv0_vg_vX = get_Multiarray_Operator(s_e->cv0_vg_vc[curved],op_inds); break;
		default: EXIT_ERROR("Unsupported: %c",node_kind); break;
	}
	const int n_vs = (int)cv0_vg_vX->op_std->ext_0;

	const struct const_Multiarray_R*const geom_coef = ( !using_restart ? s_vol->geom_coef : s_vol->geom_coef_p1 );
	struct Multiarray_R* xyz_v = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_vs,d}); // returned
	mm_NN1C_Operator_Multiarray_R(cv0_vg_vX,geom_coef,xyz_v,op_format,geom_coef->order,NULL,NULL);

	return (const struct const_Multiarray_R*) xyz_v;
}

static const struct const_Multiarray_R* constructor_xyz_f
	(const struct Simulation* sim, struct Solver_Face_T* s_face, const char node_kind, const bool using_restart)
{
	UNUSED(sim);

	struct Face* face  = (struct Face*) s_face;
	struct Volume* vol = (struct Volume*) face->neigh_info[0].volume;

	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int ind_lf   = face->neigh_info[0].ind_lf,
	          ind_href = face->neigh_info[0].ind_href,
	          curved = ( !using_restart ? vol->curved : false );
	const int p_v = ( !curved ? 1 : ((struct Solver_Volume_T*)vol)->p_ref ),
	          p_f = s_face->p_ref;

	const struct Operator* cv0_vg_fX = NULL;
	if (node_kind == 'f')
		cv0_vg_fX = get_Multiarray_Operator(s_e->cv0_vg_ff[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
	else if (node_kind == 'c')
		cv0_vg_fX = get_Multiarray_Operator(s_e->cv0_vg_fc[curved],(ptrdiff_t[]){ind_lf,ind_href,0,p_f,p_v});
	else if (node_kind == 't')
		EXIT_ADD_SUPPORT;
	else
		EXIT_UNSUPPORTED;

	const int d    = ((struct const_Element*)s_e)->d,
	          n_fs = (int)cv0_vg_fX->op_std->ext_0;

	struct Multiarray_R* xyz_f = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_fs,d}); // returned

	const char op_format = 'd';

	const struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) vol;
	const struct const_Multiarray_R*const geom_coef = ( !using_restart ? s_vol->geom_coef : s_vol->geom_coef_p1 );
	mm_NN1C_Operator_Multiarray_R(cv0_vg_fX,geom_coef,xyz_f,op_format,geom_coef->order,NULL,NULL);

	return (const struct const_Multiarray_R*) xyz_f;
}

static void compute_coef_from_val_vs
	(const struct Solver_Volume_T* s_vol, const struct const_Multiarray_T* sol_val, struct Multiarray_T* sol_coef)
{
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int p_ref = s_vol->p_ref;

	const struct Operator* vc0_vs_vs = get_Multiarray_Operator(s_e->vc0_vs_vs,(ptrdiff_t[]){0,0,p_ref,p_ref});

	resize_Multiarray_T(sol_coef,sol_val->order,sol_val->extents);
	mm_NN1C_Operator_Multiarray_T(vc0_vs_vs,sol_val,sol_coef,op_format,sol_coef->order,NULL,NULL);
}

static void compute_coef_from_val_vg
	(const struct Solver_Volume_T* s_vol, const struct const_Multiarray_T* grad_val, struct Multiarray_T* grad_coef)
{
	const char op_format = 'd';

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int p_ref = s_vol->p_ref;

	const struct Operator* vc0_vr_vr = get_Multiarray_Operator(s_e->vc0_vr_vr,(ptrdiff_t[]){0,0,p_ref,p_ref});

	resize_Multiarray_T(grad_coef,grad_val->order,grad_val->extents);
	mm_NN1C_Operator_Multiarray_T(vc0_vr_vr,grad_val,grad_coef,op_format,grad_coef->order,NULL,NULL);
}

static const struct const_Multiarray_T* constructor_nf
	(struct Solver_Face_T* s_face, const char node_kind, struct Flux_Input_T*const flux_i,
	 const struct Simulation*const sim)
{
	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;

	struct Solution_Container_T sol_cont =
		{ .ce_type = 'f', .cv_type = 'v', .node_kind = node_kind, .volume = NULL, .face = NULL, .sol = NULL, };
	struct Multiarray_T* sol_fs = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){0,0}); // destructed

	sol_cont.face = s_face;
	sol_cont.sol  = sol_fs;
	test_case->set_sol_start(sim,sol_cont);

	flux_i->s = (struct const_Multiarray_T*)sol_fs;
	flux_i->xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed

	struct Flux_T* flux = constructor_Flux_T(flux_i);
	destructor_const_Multiarray_R(flux_i->xyz);

	const struct const_Multiarray_R* normals = NULL;
	bool destruct_normals = false;
	if (node_kind == 'f') {
		normals = constructor_normals_ff(s_face); // destructed
		destruct_normals = true;
	} else if (node_kind == 'c') {
		normals = s_face->normals_fc;
	} else {
		EXIT_UNSUPPORTED;
	}

	const struct const_Multiarray_T* nf = constructor_n_dot_f(flux->f,normals); // destructed
	destructor_Flux_T(flux);

	if (destruct_normals)
		destructor_const_Multiarray_R(normals);

	destructor_Multiarray_T(sol_fs);

	return nf;
}

static void compute_coef_from_val_ff
	(const struct Solver_Face_T* s_face, const struct const_Multiarray_T* f_val, struct Multiarray_T* f_coef)
{
	const struct Face* face            = (struct Face*) s_face;
	const struct Solution_Element* s_e = &((struct Solver_Element*)face->neigh_info[0].volume->element)->s_e;

	const int ind_e = get_face_element_index(face),
	          p     = s_face->p_ref;
	const struct Operator* vc0_ff_ff = get_Multiarray_Operator(s_e->vc0_ff_ff,(ptrdiff_t[]){ind_e,ind_e,0,0,p,p});

	resize_Multiarray_T(f_coef,f_val->order,f_val->extents);
	mm_NN1C_Operator_Multiarray_T(vc0_ff_ff,f_val,f_coef,'d',f_coef->order,NULL,NULL);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for the multiarray of metric terms at the face flux nodes.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_metrics_ff
	(const struct Solver_Face_T* s_face ///< \ref Solver_Face_T.
	);

static const struct const_Multiarray_R* constructor_normals_ff (const struct Solver_Face_T* s_face)
{
	const struct const_Multiarray_R* metrics_ff = constructor_metrics_ff(s_face); // destructed

	const struct Face* face       = (struct Face*) s_face;
	const struct Volume* vol      = face->neigh_info[0].volume;
	const struct const_Element* e = (struct const_Element*) vol->element;

	struct Multiarray_R* normals_ff = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){0,0}); // returned
	compute_unit_normals(face->neigh_info[0].ind_lf,e->normals,metrics_ff,normals_ff);
	destructor_const_Multiarray_R(metrics_ff);

	transpose_Multiarray_R(normals_ff,true); // Convert to column major.

	return (struct const_Multiarray_R*) normals_ff;
}

static const struct const_Multiarray_T* constructor_n_dot_f
	(const struct const_Multiarray_T* flux, const struct const_Multiarray_R* normals)
{
	assert(flux->layout == 'C');
	assert(flux->extents[0] == normals->extents[0]);
	assert(flux->extents[1] == normals->extents[1]);

	const bool need_transpose = ( flux->layout != normals->layout);
	if (need_transpose)
		transpose_Multiarray_R((struct Multiarray_R*)normals,true);

	const int n_n  = (int)flux->extents[0],
	          d    = (int)flux->extents[1],
	          n_eq = (int)flux->extents[2];

	struct Multiarray_T* nf = constructor_zero_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_eq}); // returned

	for (int eq = 0; eq < n_eq; ++eq) {
		Type*const data_nf = get_col_Multiarray_T(eq,nf);
		for (int dim = 0; dim < d; ++dim) {
			const ptrdiff_t ind_f = compute_index_sub_container(flux->order,1,flux->extents,(ptrdiff_t[]){dim,eq});
			const Type*const data_f = &flux->data[ind_f];
			const Real*const data_n = get_col_const_Multiarray_R(dim,normals);
			for (int n = 0; n < n_n; ++n)
				data_nf[n] += data_f[n]*data_n[n];
		}
	}

	if (need_transpose)
		transpose_Multiarray_R((struct Multiarray_R*)normals,true);

	return (struct const_Multiarray_T*) nf;
}

// Level 2 ********************************************************************************************************** //

static const struct const_Multiarray_R* constructor_metrics_ff (const struct Solver_Face_T* s_face)
{
	const struct Face* face           = (struct Face*) s_face;
	const struct Volume* vol          = face->neigh_info[0].volume;
	const struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) vol;

	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int ind_lf = face->neigh_info[0].ind_lf,
	          cv     = vol->curved,
	          pv     = ( cv ? s_vol->p_ref : 1),
	          pf     = s_face->p_ref;
	const struct Operator* vv0_vm_ff = get_Multiarray_Operator(s_e->vv0_vm_ff[cv],(ptrdiff_t[]){ind_lf,0,0,pf,pv});

	const struct const_Multiarray_R* m_vm = s_vol->metrics_vm;
	return constructor_mm_NN1_Operator_const_Multiarray_R(vv0_vm_ff,m_vm,'C','d',m_vm->order,NULL);
}
