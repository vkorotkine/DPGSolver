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

#include "compute_error.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_error.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "element_solution.h"
#include "element_solver.h"
#include "face_solver.h"
#include "volume_solver.h"

#include "const_cast.h"
#include "file_processing.h"
#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Destructor for a \ref Error_CE container.
static void destructor_Error_CE
	(struct Error_CE* error_ce ///< Standard.
	);

/** \brief Compute the volume of the domain.
 *  \return See brief. */
static double compute_domain_volume
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Increment the global squared \f$L^2\f$ errors with the contribution from the current volume.
static void increment_vol_errors_l2_2
	(struct Vector_d* errors_l2_2,           ///< Holds the global squared l2 errors.
	 const struct const_Multiarray_d* err_v, ///< Holds the error values for the current volume.
	 const struct Solver_Volume* s_vol       ///< Current \ref Solver_Volume_T.
	);

/// \brief Increment the global errors with the integrated contribution from the current face.
static void increment_face_errors_integrated
	(struct Vector_d*const errors_int,            ///< Holds the integrated errors over the face.
	 const struct const_Multiarray_d*const err_f, ///< Holds the error values for the current face.
	 const struct Solver_Face*const s_face,       ///< Current \ref Solver_Face_T.
	 const int pow_exponent                       ///< The exponent to which to raise the error term.
	);

/// \brief Output the errors to the 's'erial/'p'arallel file.
static void output_errors_sp
	(const char sp_type,              ///< 's'erial/'p'arallel specifier.
	 const int error_type,            ///< The number representing the type of error being output.
	 const struct Error_CE* error_ce, ///< \ref Error_CE.
	 const struct Simulation* sim     ///< \ref Simulation.
	);

/// \brief Output the combined errors from all processors.
static void output_errors_global
	(const int error_type,            ///< The number representing the type of error being output.
	 const struct Error_CE* error_ce, ///< \ref Error_CE.
	 const struct Simulation* sim     ///< \ref Simulation.
	);

/** \brief Constructor for a derived \ref Solver_Volume_T used to compute the exact solution.
 *  \return See brief. */
static struct Solver_Volume* constructor_Solver_Volume_exact ();

/// \brief Destructor for a derived \ref Solver_Volume_T used to compute the exact solution.
static void destructor_Solver_Volume_exact
	(struct Solver_Volume* s_vol_ex ///< Standard.
	);

/// \brief Set the relevant members of a duplicate \ref Solver_Volume_T used to compute the exact solution.
static void set_Solver_Volume_exact
	(struct Solver_Volume* s_vol_ex, ///< The partially duplicated \ref Solver_Volume_T.
	 struct Solver_Volume* s_vol     ///< The \ref Solver_Volume_T.
	);

// Interface functions ********************************************************************************************** //

void output_error (const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	struct Error_CE* error_ce = test_case->constructor_Error_CE(sim);

	output_errors_sp('s',ERROR_STANDARD,error_ce,sim);
	output_errors_global(ERROR_STANDARD,error_ce,sim);

	destructor_Error_CE(error_ce);
}

void output_error_functionals (const struct Simulation*const sim)
{
	assert(sim->volumes->name == IL_VOLUME_SOLVER);
	assert(sim->faces->name   == IL_FACE_SOLVER);
	assert(list_is_derived_from("solver",'e',sim));

	const struct Test_Case*const test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (!test_case->constructor_Error_CE_functionals)
		return;

	assert(test_case->has_functional);

	struct Error_CE* error_ce = test_case->constructor_Error_CE_functionals(sim);

	output_errors_sp('s',ERROR_FUNCTIONAL,error_ce,sim);
	output_errors_global(ERROR_FUNCTIONAL,error_ce,sim);

	destructor_Error_CE(error_ce);
}

struct Error_CE* constructor_Error_CE (struct Error_CE_Helper* e_ce_h, const struct Simulation* sim)
{
	// Finalize sol_err
	struct Vector_d* sol_err = e_ce_h->sol_err;
	switch (e_ce_h->error_type) {
	case ERROR_STANDARD:
		for (int i = 0; i < sol_err->ext_0; ++i)
			sol_err->data[i] = sqrt(sol_err->data[i]/(e_ce_h->domain_volume));
		break;
	case ERROR_FUNCTIONAL:
		break; // do nothing.
	default:
		EXIT_ERROR("Unsupported: %d\n",e_ce_h->error_type);
		break;
	}

	// Construct Error_CE
	struct Error_CE* error_ce = calloc(1,sizeof *error_ce); // returned

	const_cast_ptrdiff(&error_ce->dof,compute_dof_sol_1st(sim));
	const_cast_d(&error_ce->domain_volume,e_ce_h->domain_volume);


	const_cast_c1(&error_ce->header_spec,e_ce_h->header_spec);

	int exp_o = -1;
	switch (e_ce_h->error_type) {
	case ERROR_STANDARD:
		exp_o = e_ce_h->domain_order+1;
		break;
	case ERROR_FUNCTIONAL:
		exp_o = 2*e_ce_h->domain_order;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_ce_h->error_type);
		break;
	}
	error_ce->sol_err = constructor_copy_const_Vector_d((struct const_Vector_d*)e_ce_h->sol_err); // destructed

	struct Vector_i* expected_order = constructor_empty_Vector_i(e_ce_h->n_out); // moved
	set_to_value_Vector_i(expected_order,exp_o);
	error_ce->expected_order = (const struct const_Vector_i*) expected_order; // destructed

	return error_ce;
}

struct Error_CE_Helper* constructor_Error_CE_Helper (const struct Simulation* sim, const int n_out)
{
	struct Error_CE_Helper* e_ce_h = malloc(sizeof * e_ce_h); // destructed

	const_cast_i(&e_ce_h->n_out,n_out);
	e_ce_h->domain_order = -1;
	e_ce_h->domain_volume = compute_domain_volume(sim);

	e_ce_h->header_spec = NULL;

	e_ce_h->sol_err = constructor_empty_Vector_d(n_out); // destructed
	set_to_value_Vector_d(e_ce_h->sol_err,0.0);

	e_ce_h->sol_cont = calloc(1,sizeof *e_ce_h->sol_cont); // free
	const_cast_c(&e_ce_h->sol_cont->ce_type,'v');
	const_cast_c(&e_ce_h->sol_cont->cv_type,'v');
	const_cast_c(&e_ce_h->sol_cont->node_kind,'c');
	e_ce_h->sol_cont->volume = NULL;
	e_ce_h->sol_cont->face   = NULL;
	e_ce_h->sol_cont->sol    = NULL;

	e_ce_h->s_vol[0] = NULL;
	e_ce_h->s_vol[1] = constructor_Solver_Volume_exact(); // destructed

	return e_ce_h;
}

void destructor_Error_CE_Helper (struct Error_CE_Helper* e_ce_h)
{
	destructor_Vector_d(e_ce_h->sol_err);

	free(e_ce_h->sol_cont);
	destructor_Solver_Volume_exact(e_ce_h->s_vol[1]);
	free(e_ce_h);
}

struct Error_CE_Data* constructor_Error_CE_Data
	(struct Error_CE_Helper* e_ce_h, const struct Simulation* sim)
{
	struct Error_CE_Data* e_ce_d = calloc(1,sizeof *e_ce_d); // destructed

	set_Solver_Volume_exact(e_ce_h->s_vol[1],e_ce_h->s_vol[0]);

	e_ce_d->sol[0] = constructor_sol_v(sim,e_ce_h->s_vol[0],e_ce_h->sol_cont->node_kind); // destructed
	e_ce_d->sol[1] = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0});            // destructed

	e_ce_h->sol_cont->sol    = e_ce_d->sol[1];
	e_ce_h->sol_cont->volume = e_ce_h->s_vol[1];
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	if (!using_restart())
		test_case->set_sol(sim,*(e_ce_h->sol_cont));
	else
		test_case->set_sol_start(sim,*(e_ce_h->sol_cont)); // Compare with restart file as exact solution.

	if (test_case->copy_initial_rhs) {
		e_ce_d->rhs[0] = constructor_rhs_v(sim,e_ce_h->s_vol[0],e_ce_h->sol_cont->node_kind); // destructed
		e_ce_d->rhs[1] = constructor_zero_Multiarray_d('C',2,e_ce_d->rhs[0]->extents);        // destructed
	}

	return e_ce_d;
}

void destructor_Error_CE_Data (struct Error_CE_Data* e_ce_d)
{
	for (int i = 0; i < 2; ++i) {
		destructor_Multiarray_d(e_ce_d->sol[i]);
		if (e_ce_d->rhs[i])
			destructor_Multiarray_d(e_ce_d->rhs[i]);
	}
	free(e_ce_d);
}

void increment_sol_L2 (struct Error_CE_Helper* e_ce_h, struct Error_CE_Data* e_ce_d)
{
	subtract_in_place_Multiarray_d(e_ce_d->sol[0],(const struct const_Multiarray_d*)e_ce_d->sol[1]);
	increment_vol_errors_l2_2(e_ce_h->sol_err,(const struct const_Multiarray_d*)e_ce_d->sol[0],e_ce_h->s_vol[0]);
}

void increment_sol_integrated_face (struct Error_CE_Helper*const e_ce_h, struct Error_CE_Data*const e_ce_d)
{
	subtract_in_place_Multiarray_d(e_ce_d->sol[0],(struct const_Multiarray_d*)e_ce_d->sol[1]);
	increment_face_errors_integrated(e_ce_h->sol_err,(struct const_Multiarray_d*)e_ce_d->sol[0],e_ce_h->s_face,1);
}

void increment_sol_face_L2 (struct Error_CE_Helper*const e_ce_h, struct Error_CE_Data*const e_ce_d)
{
	subtract_in_place_Multiarray_d(e_ce_d->sol[0],(struct const_Multiarray_d*)e_ce_d->sol[1]);
	increment_face_errors_integrated(e_ce_h->sol_err,(struct const_Multiarray_d*)e_ce_d->sol[0],e_ce_h->s_face,2);
}

void update_domain_order (struct Error_CE_Helper* e_ce_h)
{
	if (e_ce_h->domain_order == -1)
		e_ce_h->domain_order = e_ce_h->s_vol[0]->p_ref;
	else
		assert(e_ce_h->domain_order == e_ce_h->s_vol[0]->p_ref);
}

const char* compute_error_file_name (const int error_type, const struct Simulation*const sim)
{
	static const char* name_part = "../output/errors/";

	char l2_spec[STRLEN_MIN];
	switch (error_type) {
	case ERROR_STANDARD:
		strcpy(l2_spec,"");
		break;
	case ERROR_FUNCTIONAL:
		strcpy(l2_spec,"functional_");
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",error_type);
		break;
	}

	static char output_name[5*STRLEN_MAX];
	sprintf(output_name,"%s%s%c%s%c%s%s%s",
	        name_part,sim->pde_name,'/',sim->pde_spec,'/',
	        l2_spec,"l2_errors__",extract_name(sim->ctrl_name_full,true));
	correct_file_name_ml_p(sim->ml_p_curr[0],sim->ml_p_curr[1],output_name);
	return output_name;
}

void correct_file_name_ml_p (const int ml, const int p, char*const file_name)
{
	int index = (int)strlen(file_name);
	if (ml >= 0 && strstr(file_name,"__ml") == NULL) {
		assert(ml < 10);
		index += sprintf(file_name+index,"%s%d","__ml",ml);
	}
	if (p >= 0 && strstr(file_name,"__p") == NULL) {
		assert(p < 10);
		index += sprintf(file_name+index,"%s%d","__p",p);
	}
}

void add_rhs_Error_CE_Data (struct Error_CE_Data*const e_ce_d, const struct Simulation*const sim)
{
	struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
	assert(test_case->copy_initial_rhs == true);

	const int n_var = test_case->n_var;

	const ptrdiff_t ext_0 = e_ce_d->sol[0]->extents[0];
	for (int i = 0; i < 2; ++i) {
		struct Multiarray_d*const sol = e_ce_d->sol[i];

		const ptrdiff_t ext_1_old = sol->extents[1];
		const ptrdiff_t ext_1_new = ext_1_old+n_var;
		resize_Multiarray_d(sol,sol->order,(ptrdiff_t[]){ext_0,ext_1_new});

		const ptrdiff_t exts_rhs[] = { sol->extents[0], n_var, };
		double*const rhs_ptr = get_col_Multiarray_d(ext_1_old,sol);
		struct Multiarray_d*const sol_rhs =
			constructor_move_Multiarray_d_d('C',sol->order,exts_rhs,true,rhs_ptr); // destructed

		copy_into_Multiarray_d(sol_rhs,(struct const_Multiarray_d*)e_ce_d->rhs[i]);
		sol_rhs->owns_data = false;
		destructor_Multiarray_d(sol_rhs);

		sol->extents[1] = ext_1_new;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

///\{ \name Summation methods for combining error contributions across the domain.
#define SUMMATION_L2_VOL   201 ///< See \ref output_errors_global.
#define SUMMATION_STANDARD 202 ///< See \ref output_errors_global.
///\}

/** \brief Compute the volume of the input \ref Solver_Volume_T.
 *  \return See brief. */
static double compute_volume
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume_T.
	);

/** \brief Constructor for a \ref const_Vector_T\* holding the cubature weights multiplied by the Jacobian determinants.
 *  \return See brief. */
static const struct const_Vector_d* constructor_w_detJ
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume_T.
	);

/** \brief Constructor for a \ref const_Vector_T\* holding the cubature weights multiplied by the Jacobian determinants
 *         for the input face.
 *  \return See brief. */
static const struct const_Vector_d* constructor_w_detJ_face
	(const struct Solver_Face*const s_face ///< \ref Solver_Face_T.
	);

static void destructor_Error_CE (struct Error_CE* error_ce)
{
	destructor_const_Vector_d(error_ce->sol_err);
	destructor_const_Vector_i(error_ce->expected_order);
	free(error_ce);
}

static double compute_domain_volume (const struct Simulation* sim)
{
	double domain_volume = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
//if(!((struct Volume*)curr)->boundary)
//	continue;
		domain_volume += compute_volume((struct Solver_Volume*) curr);

	}
	return domain_volume;
}

static void increment_vol_errors_l2_2
	(struct Vector_d* errors_l2_2, const struct const_Multiarray_d* err_v, const struct Solver_Volume* s_vol)
{
	assert(errors_l2_2->ext_0 == err_v->extents[1]);

	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	const ptrdiff_t n_out = errors_l2_2->ext_0;
	for (int i = 0; i < n_out; ++i) {
		const double* err_data = get_col_const_Multiarray_d(i,err_v);
		for (int j = 0; j < ext_0; ++j)
			errors_l2_2->data[i] += w_detJ->data[j]*err_data[j]*err_data[j];
	}
	destructor_const_Vector_d(w_detJ);
}

static void increment_face_errors_integrated
	(struct Vector_d*const errors_int, const struct const_Multiarray_d*const err_f,
	 const struct Solver_Face*const s_face, const int pow_exponent)
{
	assert(errors_int->ext_0 == err_f->extents[1]);
	assert(pow_exponent == 1 || pow_exponent == 2); // Think before adding further support.

	const struct const_Vector_d*const w_detJ = constructor_w_detJ_face(s_face); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	const ptrdiff_t n_out = errors_int->ext_0;
	for (int i = 0; i < n_out; ++i) {
		const double* err_data = get_col_const_Multiarray_d(i,err_f);
		for (int j = 0; j < ext_0; ++j)
			errors_int->data[i] += w_detJ->data[j]*pow(err_data[j],pow_exponent);
	}
	destructor_const_Vector_d(w_detJ);
}

static void output_errors_sp
	(const char sp_type, const int error_type, const struct Error_CE* error_ce, const struct Simulation* sim)
{
	const ptrdiff_t n_out         = error_ce->sol_err->ext_0;
	const double*const error_data = error_ce->sol_err->data;

	const char* output_name = compute_error_file_name(error_type,sim);
	FILE* sp_file = fopen_sp_output_file(sp_type,output_name,"txt",sim->mpi_rank); // closed

	if (sp_type == 'p') {
		fprintf(sp_file,"n_out: %td\n",n_out);
		fprintf(sp_file,"expected_orders: ");

		for (int i = 0; i < n_out; ++i)
			fprintf(sp_file,"%2d ",error_ce->expected_order->data[i]);
		fprintf(sp_file,"\n\n");

		fprintf(sp_file,"%-14s%s\n","h",error_ce->header_spec);
		fprintf(sp_file,"%-14.4e",1.0/pow((double)error_ce->dof,1.0/DIM));
	} else if (sp_type == 's') {
		fprintf(sp_file,"%-10s%-14s%s\n","dof","vol",error_ce->header_spec);
		fprintf(sp_file,"%-10td%-14.4e",error_ce->dof,error_ce->domain_volume);
	}

	for (int i = 0; i < n_out; ++i)
		fprintf(sp_file,"%-14.4e",error_data[i]);

	fclose(sp_file);
}

static void output_errors_global (const int error_type, const struct Error_CE* error_ce, const struct Simulation* sim)
{
	if (!(sim->mpi_rank == 0))
		return;
	MPI_Barrier(MPI_COMM_WORLD);

	int summation_method = -1;
	switch (error_type) {
	case ERROR_STANDARD:
		summation_method = SUMMATION_L2_VOL;
		break;
	case ERROR_FUNCTIONAL:
		summation_method = SUMMATION_STANDARD;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",error_type);
		break;
	}

	const char* file_name = compute_error_file_name(error_type,sim);

	const ptrdiff_t n_out = error_ce->sol_err->ext_0;

	ptrdiff_t dof_g = 0;
	double domain_volume_g = 0.0;
	struct Vector_d* sol_error_g = constructor_empty_Vector_d(n_out); // moved
	set_to_value_Vector_d(sol_error_g,0.0);

	for (int mpi_rank = 0; mpi_rank < sim->mpi_size; ++mpi_rank) {
		FILE* s_file = fopen_sp_input_file('s',file_name,"txt",mpi_rank);

		skip_lines(s_file,1);
		char line[STRLEN_MAX];
		char* line_ptr[1] = {line};
		fgets_checked(line,sizeof(line),s_file);

		int dof = 0;
		read_line_values_i(line_ptr,1,&dof,false);
		dof_g += dof;

		double domain_volume = 0.0;
		read_skip_d_1(line,1,&domain_volume,1);
		domain_volume_g += domain_volume;

		double sol_e[n_out];
		read_skip_d_1(line,2,sol_e,(int)n_out);
		for (int i = 0; i < sol_error_g->ext_0; ++i) {
			switch (summation_method) {
				case SUMMATION_L2_VOL:   sol_error_g->data[i] += sol_e[i]*sol_e[i]*domain_volume; break;
				case SUMMATION_STANDARD: sol_error_g->data[i] += sol_e[i];                        break;
				default: EXIT_ERROR("Unsupported: %d\n",summation_method);                        break;
			}
		}
		fclose(s_file);
	}

	for (int i = 0; i < sol_error_g->ext_0; ++i) {
		switch (summation_method) {
			case SUMMATION_L2_VOL:   sol_error_g->data[i] = sqrt(sol_error_g->data[i]/domain_volume_g); break;
			case SUMMATION_STANDARD: sol_error_g->data[i] = fabs(sol_error_g->data[i]);                 break;
			default: EXIT_ERROR("Unsupported: %d\n",summation_method);                                  break;
		}
	}

	struct Error_CE* error_ce_g = calloc(1,sizeof *error_ce_g); // destructed

	const_cast_d(&error_ce_g->domain_volume,domain_volume_g);
	error_ce_g->sol_err = (const struct const_Vector_d*) sol_error_g; // keep
	error_ce_g->expected_order = constructor_copy_const_Vector_i_i(n_out,error_ce->expected_order->data); // keep
	const_cast_c1(&error_ce_g->header_spec,error_ce->header_spec);
	const_cast_ptrdiff(&error_ce_g->dof,dof_g);

	output_errors_sp('p',error_type,error_ce_g,sim);

	destructor_Error_CE(error_ce_g);
}

static struct Solver_Volume* constructor_Solver_Volume_exact ()
{
	struct Solver_Volume* s_vol_ex = calloc(1,sizeof *s_vol_ex); // returned

	s_vol_ex->sol_coef = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){0,0}); // destructed

	return s_vol_ex;
}

static void destructor_Solver_Volume_exact (struct Solver_Volume* s_vol_ex)
{
	destructor_Multiarray_d(s_vol_ex->sol_coef);
	free(s_vol_ex);
}

static void set_Solver_Volume_exact (struct Solver_Volume* s_vol_ex, struct Solver_Volume* s_vol)
{
	struct Volume* vol_ex = (struct Volume*) s_vol_ex,
	             * vol    = (struct Volume*) s_vol;

	const_cast_b(&vol_ex->curved,vol->curved);
	const_cast_const_Element(&vol_ex->element,vol->element);

	const_cast_i(&s_vol_ex->p_ref,s_vol->p_ref);
	const_constructor_move_const_Multiarray_d(&s_vol_ex->geom_coef,s_vol->geom_coef);
	const_constructor_move_const_Multiarray_d(&s_vol_ex->geom_coef_p1,s_vol->geom_coef_p1);
}

// Level 1 ********************************************************************************************************** //

static double compute_volume (const struct Solver_Volume* s_vol)
{
	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	double volume = 0.0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_detJ->data[i];

	destructor_const_Vector_d(w_detJ);

	return volume;
}

static const struct const_Vector_d* constructor_w_detJ (const struct Solver_Volume* s_vol)
{
	struct Volume* vol = (struct Volume*)s_vol;
	const struct Solution_Element* s_e = &((struct Solver_Element*)vol->element)->s_e;

	const int curved = vol->curved,
	          p      = s_vol->p_ref;
	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(s_e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	assert(w_vc->ext_0 == jacobian_det_vc.ext_0);

	const ptrdiff_t ext_0 = w_vc->ext_0;

	struct Vector_d* w_detJ = constructor_empty_Vector_d(ext_0); // returned

	for (int i = 0; i < ext_0; ++i)
		w_detJ->data[i] = w_vc->data[i]*jacobian_det_vc.data[i];

	return (const struct const_Vector_d*) w_detJ;
}

static const struct const_Vector_d* constructor_w_detJ_face (const struct Solver_Face*const s_face)
{
	const struct const_Vector_d*const w_fc = get_operator__w_fc__s_e(s_face);
	const struct const_Vector_d jacobian_det_fc = interpret_const_Multiarray_as_Vector_d(s_face->jacobian_det_fc);
	assert(w_fc->ext_0 == jacobian_det_fc.ext_0);

	const ptrdiff_t ext_0 = w_fc->ext_0;

	struct Vector_d* w_detJ = constructor_empty_Vector_d(ext_0); // returned
	for (int i = 0; i < ext_0; ++i)
		w_detJ->data[i] = w_fc->data[i]*jacobian_det_fc.data[i];

	return (const struct const_Vector_d*) w_detJ;
}
