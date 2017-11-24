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
#include <string.h>
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
#include "core.h"
#include "file_processing.h"
#include "geometry.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

///\{ \name Flag for whether the convergence orders should be displayed for these tests.
#define DISPLAY_CONV 1
///\}

///\{ \name Parameters relating to which solutions to output to paraview for visualization.
#define ORDER_VIS_CONV_P      2
#define ORDER_VIS_CONV_ML_MAX 3
///\}

/// \brief Check the convergence orders of the errors for the simulations performed.
static void check_convergence_orders
	(struct Test_Info*const test_info,                 ///< \ref Test_Info.
	 const struct Integration_Test_Info* int_test_info ///< \ref Integration_Test_Info.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the solution convergence (\ref test_integration_convergence.c).
 *  \return 0 on success (when the solution converges at the expected rate).
 */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(nargc == 3,"Invalid number of input arguments");

	const char* petsc_options_name = set_petsc_options_name(argv[2]);
	PetscInitialize(&nargc,&argv,petsc_options_name,PETSC_NULL);

	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim = NULL;
	for (int ml = ml_ref[0], ml_prev = ml-1; ml <= ml_ref[1]; ++ml) {
	for (int p = p_ref[0], p_prev = p-1; p <= p_ref[1]; ++p) {
		const int adapt_type = int_test_info->adapt_type;
		const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr); // destructed

		if ((ml == ml_ref[0]) && (p == p_ref[0]))
			const_cast_d(&int_test_info->conv_order_discount,sim->test_case->conv_order_discount);

		solve_for_solution(sim);

		if (p == ORDER_VIS_CONV_P && ml <= ORDER_VIS_CONV_ML_MAX) {
			output_visualization(sim,VIS_GEOM_EDGES);
			output_visualization(sim,VIS_SOLUTION);
		}

		output_error(sim);

		if (DISPLAY_CONV)
			printf("ml, p, dof: %d %d %td\n",ml,p,compute_dof(sim));

		p_prev  = p;
		ml_prev = ml;
		structor_simulation(&sim,'d',adapt_type,p,ml,p_prev,ml_prev,NULL);
	}}

	check_convergence_orders(&test_info,int_test_info);
	destructor_Integration_Test_Info(int_test_info);

	PetscFinalize();
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the number of errors to be read for the convergence order test.
 *  \return See brief. */
static int compute_n_err
	(const char* input_name ///< \ref fopen_sp_input_file :: name_part.
	);

/** \brief Return whether the convergence orders are in the expected range.
 *  \return See brief. */
static bool attained_expected_conv_orders
	(const double discount,                            ///< Allowable discount from the expected conv. orders.
	 const ptrdiff_t* extents_e,                       ///< The extents of the error data.
	 const struct Multiarray_d* conv_orders,           ///< The container for the conv. order data.
	 const struct Multiarray_i* exp_orders,            ///< The container of expected conv. order data.
	 const struct Integration_Test_Info* int_test_info ///< \ref Integration_Test_Info.
	);

static void check_convergence_orders
	(struct Test_Info*const test_info, const struct Integration_Test_Info* int_test_info)
{
	UNUSED(test_info);
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
		const char* input_name_curr = set_file_name_curr(ADAPT_0,p,ml,true,input_name);

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

	if (DISPLAY_CONV) {
		printf("h:\n");
		print_Multiarray_d(h);
		printf("L2 errors:\n");
		print_Multiarray_d(l2_err);
		printf("Convergence orders:\n");
		print_Multiarray_d(conv_orders);
	}

	const double disc = 0.125;
	bool pass   = attained_expected_conv_orders(disc,extents,conv_orders,ex_ord,int_test_info);
	bool pass_d = attained_expected_conv_orders(
		disc+int_test_info->conv_order_discount,extents,conv_orders,ex_ord,int_test_info);

	destructor_Multiarray_i(ex_ord);
	destructor_Multiarray_d(l2_err);
	destructor_Multiarray_d(h);
	destructor_Multiarray_d(conv_orders);

	if (pass_d && !pass) {
		test_print_warning(test_info,"Only passing with the convergence order discount.");
		pass = true;
	}
	assert_condition(pass);
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

static bool attained_expected_conv_orders
	(const double discount, const ptrdiff_t* extents_e, const struct Multiarray_d* conv_orders,
	 const struct Multiarray_i* exp_orders, const struct Integration_Test_Info* int_test_info)
{
	bool pass = true;

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;
	const int n_err   = extents_e[2];
	const ptrdiff_t base_e = extents_e[0]*extents_e[1];

	const int ml = ml_ref[1];
	for (int p = p_ref[0]; p <= p_ref[1]; ++p) {
	for (int n = 0; n < n_err; ++n) {
		const ptrdiff_t ind_e = n*base_e + p*extents_e[0] + ml;
		if (isnan(conv_orders->data[ind_e]) || (conv_orders->data[ind_e] < (exp_orders->data[ind_e]-discount))) {
			pass = false;
			break;
		}
	}}
	return pass;
}
