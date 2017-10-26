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

#include "test_integration_convergence.h"

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_intrusive.h"
#include "definitions_visualization.h"

#include "test_base.h"
#include "test_integration.h"

#include "computational_elements.h"

#include "multiarray.h"

#include "compute_error.h"
#include "const_cast.h"
#include "file_processing.h"
#include "geometry.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

///\{ \name Parameters relating to which solutions to output to paraview for visualization.
#define ORDER_VIS_CONV_P      3
#define ORDER_VIS_CONV_ML_MAX 3
///\}

/// \brief 'Con'/'De'strucor for a \ref Simulation depending on `adapt_type`.
static void structor_simulation
	(struct Simulation** sim,   ///< Pointer to the \ref Simulation.
	 const char mode,           ///< Mode of operation. Options: 'c'onstructor, 'd'estructor.
	 const int adapt_type,      ///< \ref Integration_Test_Info::adapt_type.
	 const int p,               ///< The order of the current simulation.
	 const int ml,              ///< The mesh level of the current simulation.
	 const int p_prev,          ///< The order of the previous simulation.
	 const int ml_prev,         ///< The mesh level of the previous simulation.
	 const char*const ctrl_name ///< The name of the control file (used for \ref constructor_Simulation).
	);

/** \brief Set the name of the current file to be used based on the adaptation type, order and mesh level.
 *  \return See brief. */
static const char* set_file_name_curr
	(const int adapt_type,      ///< \ref Integration_Test_Info::adapt_type.
	 const int p,               ///< The order of the current simulation.
	 const int ml,              ///< The mesh level of the current simulation.
	 const char*const file_name ///< The name of the file.
	);

/// \brief Check the convergence orders of the errors for the simulations performed.
static void check_convergence_orders
	(struct Test_Info*const test_info,                 ///< \ref Test_Info.
	 const struct Integration_Test_Info* int_test_info ///< \ref Integration_Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_integration_convergence (struct Test_Info*const test_info, const char*const ctrl_name)
{
	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	for (int p = p_ref[0], p_prev = p; p <= p_ref[1]; ++p) {
	for (int ml = ml_ref[0], ml_prev = ml; ml <= ml_ref[1]; ++ml) {
		const int adapt_type = int_test_info->adapt_type;
		const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,ctrl_name);
		struct Simulation* sim = NULL;
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr);

		constructor_derived_computational_elements(sim,IL_SOLVER); // destructed
		solve_for_solution(sim);

		if (p == ORDER_VIS_CONV_P && ml <= ORDER_VIS_CONV_ML_MAX) {
			output_visualization(sim,VIS_GEOM_EDGES);
			output_visualization(sim,VIS_SOLUTION);
		}

		output_error(sim);
		destructor_derived_computational_elements(sim,IL_BASE);

		const_cast_b(&int_test_info->display_progress,sim->test_case->display_progress);
		if (int_test_info->display_progress)
			printf("ml, p, dof: %d %d %td\n",ml,p,compute_dof(sim));

		p_prev  = p;
		ml_prev = ml;
		structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL);
	}}

	check_convergence_orders(test_info,int_test_info);
	destructor_Integration_Test_Info(int_test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the number of errors to be read for the convergence order test.
 *  \return See brief. */
static int compute_n_err
	(const char* input_name ///< \ref fopen_sp_input_file :: name_part.
	);

static void structor_simulation
	(struct Simulation** sim, const char mode, const int adapt_type, const int p, const int ml, const int p_prev,
	 const int ml_prev, const char*const ctrl_name)
{
	assert(mode == 'c' || mode == 'd');

	switch (adapt_type) {
	case ADAPT_0:
		if (mode == 'c') {
			*sim = constructor_Simulation(ctrl_name); // destructed
		} else if (mode == 'd') {
			destructor_Simulation(*sim);
		}
		break;
	case ADAPT_P:
		EXIT_ADD_SUPPORT;
		if (ml != ml_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			; // p-adapt
		}
		break;
	case ADAPT_H:
		EXIT_ADD_SUPPORT;
		if (p != p_prev) {
			structor_simulation(sim,mode,ADAPT_0,p,ml,p_prev,ml_prev,ctrl_name);
		} else {
			; // h-adapt
		}
		break;
	case ADAPT_HP:
		EXIT_ADD_SUPPORT;
		if (ml != ml_prev)
			{ ; } // h-adapt
		if (p != p_prev)
			{ ; } // p-adapt
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}
}

static const char* set_file_name_curr (const int adapt_type, const int p, const int ml, const char*const file_name)
{
	static char file_name_curr[STRLEN_MAX] = { 0, };
	strcpy(file_name_curr,file_name);

	assert(p  >= 0 && p  <= 9); // Constraints are only imposed because of the currently limited string processing.
	assert(ml >= 0 && ml <= 9);

	char* index = NULL;
	switch (adapt_type) {
	case ADAPT_0:
		index = strstr(file_name_curr,"__ml");
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		index = strstr(file_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
		break;
	case ADAPT_P:
		assert(strstr(file_name_curr,"__ml") == NULL);
		index = strstr(file_name_curr,"__p");
		index[3] = '0'+p;
		assert(!isdigit(index[4]));
		break;
	case ADAPT_H:
		assert(strstr(file_name_curr,"__p") == NULL);
		index = strstr(file_name_curr,"__ml");
		index[4] = '0'+ml;
		assert(!isdigit(index[5]));
		break;
	case ADAPT_HP:
		assert(strstr(file_name_curr,"__ml") == NULL);
		assert(strstr(file_name_curr,"__p") == NULL);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",adapt_type);
		break;
	}

	return file_name_curr;
}

static void check_convergence_orders
	(struct Test_Info*const test_info, const struct Integration_Test_Info* int_test_info)
{
	struct Simulation* sim = constructor_Simulation(int_test_info->ctrl_name); // destructed
	if (sim->mpi_rank) {
		destructor_Simulation(sim);
		return;
	}

	const char* input_name = compute_error_file_name(sim);
	destructor_Simulation(sim);

	const int n_err = compute_n_err(input_name);

	ptrdiff_t extents[]    = { int_test_info->ml[1]+1, int_test_info->p_ref[1]+1, n_err, };
	const ptrdiff_t base_e = extents[0]*extents[1];

	struct Multiarray_i* ex_ord = constructor_zero_Multiarray_i('C',3,extents);     // destructed
	struct Multiarray_d* l2_err = constructor_zero_Multiarray_d('C',3,extents);     // destructed
	struct Multiarray_d* h      = constructor_zero_Multiarray_d('C',2,&extents[0]); // destructed

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;
	for (int p = p_ref[0]; p <= p_ref[1]; ++p) {
	for (int ml = ml_ref[0]; ml <= ml_ref[1]; ++ml) {
		const char* input_name_curr = set_file_name_curr(ADAPT_0,p,ml,input_name);

		char line[STRLEN_MAX];
		FILE* p_file = fopen_sp_input_file('p',input_name_curr,"txt",0); // closed

		skip_lines(p_file,1);
		fgets(line,sizeof(line),p_file);
		int data_i[n_err];
		read_skip_i_1(line,1,data_i,n_err);

		skip_lines(p_file,2);
		fgets(line,sizeof(line),p_file);

		double data_d[n_err+1];
		read_skip_d_1(line,0,data_d,n_err+1);
		fclose(p_file);

		const ptrdiff_t ind_h = compute_index_sub_container(h->order,0,h->extents,(ptrdiff_t[]){ml,p});
		h->data[ind_h] = data_d[0];

		for (int n = 0; n < n_err; ++n) {
			const ptrdiff_t ind_e = n*base_e + p*extents[0] + ml;
			l2_err->data[ind_e] = data_d[n+1];
			ex_ord->data[ind_e] = data_i[n];
		}
	}}

	struct Multiarray_d* conv_orders = constructor_zero_Multiarray_d('C',3,extents); // destructed

	for (int p = p_ref[0]; p <= p_ref[1]; ++p) {
	for (int ml = ml_ref[0]+1; ml <= ml_ref[1]; ++ml) {
	for (int n = 0; n < n_err; ++n) {
		const ptrdiff_t ind_c_h = p*extents[0] + ml-1,
		                ind_f_h = p*extents[0] + ml,
		                ind_c_e = n*base_e + ind_c_h,
		                ind_f_e = n*base_e + ind_f_h;

		double* data_l2 = l2_err->data,
		      * data_h  = h->data;

		conv_orders->data[ind_f_e] =
			log10(data_l2[ind_f_e]/data_l2[ind_c_e])/log10(data_h[ind_f_h]/data_h[ind_c_h]);
	}}}

	if (int_test_info->display_progress) {
		printf("h:\n");
		print_Multiarray_d(h);
		printf("L2 errors:\n");
		print_Multiarray_d(l2_err);
		printf("Convergence orders:\n");
		print_Multiarray_d(conv_orders);
	}

	bool pass        = true;
	const double tol = 0.125;

	const int ml = ml_ref[1];
	for (int p = p_ref[0]; p <= p_ref[1]; ++p) {
	for (int n = 0; n < n_err; ++n) {
		const ptrdiff_t ind_e = n*base_e + p*extents[0] + ml;
		if (isnan(conv_orders->data[ind_e]) || (conv_orders->data[ind_e] < (ex_ord->data[ind_e]-tol))) {
			pass = false;
			break;
		}
	}}

	destructor_Multiarray_i(ex_ord);
	destructor_Multiarray_d(l2_err);
	destructor_Multiarray_d(h);

	destructor_Multiarray_d(conv_orders);

	char test_name[STRLEN_MAX];
	sprintf(test_name,"%s%s","Conv Orders - ",int_test_info->ctrl_name);
	test_increment_and_print_name(test_info,pass,test_name);
}

// Level 1 ********************************************************************************************************** //

static int compute_n_err (const char* input_name)
{
	FILE* p_file = fopen_sp_input_file('p',input_name,"txt",0); // closed

	char line[STRLEN_MAX];
	fgets(line,sizeof(line),p_file);

	int n_err = 0;
	read_skip_i(line,&n_err);
	fclose(p_file);

	return n_err;
}
