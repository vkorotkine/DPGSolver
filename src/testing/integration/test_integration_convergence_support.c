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

#include "test_integration_convergence_support.h"

#include <assert.h>
#include <string.h>
#include <math.h>

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_error.h"
#include "definitions_visualization.h"

#include "test_base.h"
#include "test_integration.h"

#include "matrix.h"
#include "multiarray.h"

#include "const_cast.h"
#include "core.h"
#include "file_processing.h"
#include "restart_writers.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

/// The discount in the convergence order value for which the convergence may still be considered optimal (p+1).
#define ACCEPTABLE_DISCOUNT 0.165

#define DISPLAY_CONV 1 ///< Flag for whether the convergence orders should be displayed for these tests.

///\{ \name Parameters relating to which solutions to output to paraview for visualization.
#define ORDER_VIS_CONV_P      2
#define ORDER_VIS_CONV_ML_MAX 4
#define DISPLAY_GEOM          0 ///< Flag for whether the geometry should be output.
///\}

///\{ \name Parameters relating to maximum allowable mesh level and order for convergence order testing.
#define ML_MAX 4
#define P_MAX  2
///\}

/** \brief Set \ref Integration_Test_Info::conv_order_discount to the value specified for the test
 *         case (from the input file) or to 0.0 otherwise. */
static void set_convergence_order_discount
	(struct Integration_Test_Info*const int_test_info ///< \ref Integration_Test_Info.
	);

/// \brief Check the convergence orders of the errors for the simulations performed.
static void check_convergence_orders
	(const int error_type,                                   ///< The type of error.
	 bool*const pass,                                        ///< To be set based on the result of the test.
	 struct Test_Info*const test_info,                       ///< \ref Test_Info.
	 const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const struct Simulation*const sim                       ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void run_convergence_order_study (int argc, char** argv, const int conv_study_type)
{
	const char* ctrl_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(ctrl_name);
	if (argc == 4)
		int_test_info->conv_study_extension = argv[3];

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim = NULL;
	const char type_rc = 'r';

	int ml_prev = ml_ref[0]-1,
	    p_prev  = p_ref[0]-1;

	bool ignore_static = false;
	int ml_max = ml_ref[1];
	switch (conv_study_type) {
	case CONV_STUDY_SOLVE:
		break; // Do nothing
	case CONV_STUDY_SOLVE_NO_CHECK: // fallthrough
	case CONV_STUDY_RESTART:
		ignore_static = true;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",conv_study_type);
		break;
	}

	for (int ml = ml_ref[0]; ml <= ml_max; ++ml) {
	for (int p = p_ref[0]; p <= p_ref[1]; ++p) {
		const int adapt_type = int_test_info->adapt_type;
		const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);
		structor_simulation(&sim,'c',adapt_type,p,ml,p_prev,ml_prev,ctrl_name_curr,type_rc,ignore_static); // d.
		ignore_static = false;

		switch (conv_study_type) {
		case CONV_STUDY_SOLVE:          // fallthrough
		case CONV_STUDY_SOLVE_NO_CHECK:
			switch (get_set_method(NULL)) {
			case METHOD_DG: case METHOD_DPG: case METHOD_OPG: case METHOD_OPGC0:
				solve_for_solution(sim);
				break;
			case METHOD_L2_PROJ:
				set_initial_solution(sim);
				break; // do nothing.
			default:
				EXIT_ERROR("Unsupported: %d\n",get_set_method(NULL));
				break;
			}
			break;
		case CONV_STUDY_RESTART: {
			assert(using_restart() == true);
			set_initial_solution(sim);
			set_to_zero_residual(sim);
			struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
			test_case->constructor_Error_CE = test_case->constructor_Error_CE_restart_test;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",conv_study_type);
			break;
		}

		if (p == ORDER_VIS_CONV_P && ml <= ORDER_VIS_CONV_ML_MAX) {
			output_visualization(sim,VIS_GEOM_EDGES);
			output_visualization(sim,VIS_SOLUTION);
			output_visualization(sim,VIS_GEOM_VOLUMES);
			output_visualization(sim,VIS_NORMALS);
		}

		output_error(sim);
		output_error_functionals(sim);

		if (DISPLAY_CONV)
			printf("\ntest_integration_convergence (ml, p, dof): %d %d %td\n\n\n",ml,p,compute_dof(sim));

		if ((ml == ml_max) && (p == p_ref[1])) {
			output_restart(sim);

			set_convergence_order_discount(int_test_info);
			bool pass = true;
			switch (conv_study_type) {
			case CONV_STUDY_SOLVE:   // fallthrough
			case CONV_STUDY_RESTART:
				check_convergence_orders(ERROR_STANDARD,&pass,&test_info,int_test_info,sim);
				check_convergence_orders(ERROR_FUNCTIONAL,&pass,&test_info,int_test_info,sim);
				break;
			case CONV_STUDY_SOLVE_NO_CHECK:
				break; // do nothing.
			default:
				EXIT_ERROR("Unsupported: %d\n",conv_study_type);
				break;
			}
			assert_condition(pass);

			structor_simulation(&sim,'d',ADAPT_0,p,ml,p_prev,ml_prev,NULL,type_rc,ignore_static);
		}

		p_prev  = p;
		ml_prev = ml;
	}}
	destructor_Integration_Test_Info(int_test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for convergence order related data for each mesh level and polynomial order.
struct Conv_Order_Data {
	const struct const_Multiarray_d* h,           ///< The multiarray of mesh spacing.
	                               * l2_err,      ///< The multiarray of L2 errors.
	                               * conv_orders; ///< The multiarray of convergence orders.
	const struct const_Multiarray_i* cases_run,   ///< The multiarray of flags indicating which cases were run.
	                               * ex_ord;      ///< The multiarray of expected orders.

	const char*const* var_names; ///< Names of the variables for which the error and convergence orders are provided.
};

/** \brief Copy the files from the $BUILD/output/error/... subdirectory to the $BUILD/output/results/... subdirectory in
 *  the case where a convergence study is being performed. */
static void copy_error_files_for_conv_study
	(const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const int error_type,                                   ///< Defined for \ref compute_error_file_name.
	 const struct Simulation*const sim                       ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Conv_Order_Data container.
 *  \return See brief. */
static struct Conv_Order_Data* constructor_Conv_Order_Data
	(const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const int error_type,                                   ///< Current error type.
	 const struct Simulation*const sim                       ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Conv_Order_Data container.
static void destructor_Conv_Order_Data
	(const struct Conv_Order_Data*const cod ///< Standard.
	);

/** \brief Compute the number of errors to be read for the convergence order test.
 *  \return See brief. */
static int compute_n_err
	(const char* input_name ///< \ref fopen_sp_input_file :: name_part.
	);

/** \brief Return whether the convergence orders are in the expected range.
 *  \return See brief. */
static bool attained_expected_conv_orders
	(const double discount,                             ///< Allowable discount from the expected conv. orders.
	 const struct const_Multiarray_d*const conv_orders, ///< The container for the conv. order data.
	 const struct const_Multiarray_i*const exp_orders,  ///< The container of expected conv. order data.
	 const struct Integration_Test_Info* int_test_info  ///< \ref Integration_Test_Info.
	);

/// \brief Output the combined results from all runs currently stored in the \ref Conv_Order_Data container.
static void output_combined_results
	(const struct Conv_Order_Data*const cod,                 ///< \ref Conv_Order_Data.
	 const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const int error_type,                                   ///< Defined for \ref compute_error_file_name.
	 const struct Simulation*const sim                       ///< \ref Simulation.
	);

static void set_convergence_order_discount (struct Integration_Test_Info*const int_test_info)
{
	char line[STRLEN_MAX];
	bool found_discount = false;

	FILE* input_file = fopen_input('t',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"conv_order_discount")) {
			found_discount = true;
			read_skip_const_d(line,&int_test_info->conv_order_discount,1,false);
		}
	}
	fclose(input_file);

	if (!found_discount)
		const_cast_d(&int_test_info->conv_order_discount,0.0);
}

static void check_convergence_orders
	(const int error_type, bool*const pass, struct Test_Info*const test_info,
	 const struct Integration_Test_Info*const int_test_info, const struct Simulation* sim)
{
	UNUSED(test_info);
	if (sim->mpi_rank)
		return;

	switch (error_type) {
	case ERROR_STANDARD:
		break; // do nothing
	case ERROR_FUNCTIONAL: {
		const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;
		if (!test_case->has_functional)
			return;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d.\n",error_type);
		break;
	}

	copy_error_files_for_conv_study(int_test_info,error_type,sim);

	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;
	const char* input_name = compute_error_file_name(error_type,sim);

	const char*const input_name_curr = set_file_name_curr(ADAPT_0,p_ref[0],ml_ref[0],true,input_name);
	FILE* p_file = fopen_sp_input_file('p',input_name_curr,"txt",0); // closed
	if (p_file == NULL) {
		fclose(p_file);
		return;
	}

	const struct Conv_Order_Data*const cod = constructor_Conv_Order_Data(int_test_info,error_type,sim); // dest.
	const struct const_Multiarray_i*const ex_ord      = cod->ex_ord;
	const struct const_Multiarray_d*const conv_orders = cod->conv_orders;
	if (DISPLAY_CONV) {
		printf("h:\n");
		print_const_Multiarray_d(cod->h);
		printf("L2 errors:\n");
		print_const_Multiarray_d(cod->l2_err);
		printf("Convergence orders:\n");
		print_const_Multiarray_d(conv_orders);
		printf("------------------------------------------------------------------------------------------\n\n\n");
	}

	if (*pass) {
		*pass       = attained_expected_conv_orders(ACCEPTABLE_DISCOUNT,conv_orders,ex_ord,int_test_info);
		bool pass_d = attained_expected_conv_orders(
			ACCEPTABLE_DISCOUNT+int_test_info->conv_order_discount,conv_orders,ex_ord,int_test_info);

		if (pass_d && !*pass) {
			test_print_warning(test_info,"Only passing with the convergence order discount.");
			*pass = true;
		}
	}

	output_combined_results(cod,int_test_info,error_type,sim);
	destructor_Conv_Order_Data(cod);
}

// Level 1 ********************************************************************************************************** //

/** \brief Return a statically allocated pointer to the range of mesh levels or polynomial orders relating to the
 *         convergence order test.
 *  \return See brief. */
static const int* get_conv_order_range
	(const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const char mp_type                                      ///< Parameter indicator. Options: 'm'l, 'p'.
	);

/** \brief Return a statically allocated `char*` containing the root of the error input file names.
 *  \return See brief. */
static const char* compute_error_input_name_root
	(const struct Integration_Test_Info*const int_test_info, ///< \ref Integration_Test_Info.
	 const int error_type,                                   ///< Defined for \ref compute_error_file_name.
	 const bool use_default,                                 ///< Flag for whether the default should be returned.
	 const struct Simulation*const sim                       ///< Defined for \ref compute_error_file_name.
	);

static void copy_error_files_for_conv_study
	(const struct Integration_Test_Info*const int_test_info, const int error_type, const struct Simulation*const sim)
{
	if (!int_test_info->conv_study_extension)
		return;

	char command[3*STRLEN_MAX];
	const char*const input_name_i = compute_error_input_name_root(int_test_info,error_type,true,sim),
	          *const input_name_o = compute_error_input_name_root(int_test_info,error_type,false,sim);

	// Delete any existing files
	char input_name_o_root[STRLEN_MAX];
	const size_t len = strlen(input_name_o);
	size_t i = 0;
	for ( ; i < len-9; ++i)
		input_name_o_root[i] = input_name_o[i];
	input_name_o_root[i] = 0;
	sprintf(command,"%s %s%s","rm",input_name_o_root,"*");
	if (system(command))
		{ ; } // Do nothing (The directory may not exist yet).

	const int*const p_range  = get_conv_order_range(int_test_info,'p'),
	         *const ml_range = get_conv_order_range(int_test_info,'m');
	for (int p = p_range[0]; p <= p_range[1]; ++p) {
	for (int ml = ml_range[0]; ml <= ml_range[1]; ++ml) {
		char input_name_curr_i[STRLEN_MAX],
		     input_name_curr_o[STRLEN_MAX];
		strcpy(input_name_curr_i,set_file_name_curr(ADAPT_0,p,ml,true,input_name_i));
		strcpy(input_name_curr_o,set_file_name_curr(ADAPT_0,p,ml,true,input_name_o));

		FILE* p_file = fopen_sp_input_file_unchecked('p',input_name_curr_i,"txt",0); // closed
		if (!p_file)
			continue;
		fclose(p_file);

		mkdir_p_given_file_name(input_name_curr_o);
		sprintf(command,"%s %s%s %s%s","cp",input_name_curr_i,"_p.txt",input_name_curr_o,"_p.txt");
		if (system(command))
			EXIT_ERROR("Problem with system call.");
	}}
}

static int compute_n_err (const char* input_name)
{
	FILE* p_file = fopen_sp_input_file('p',input_name,"txt",0); // closed

	char line[STRLEN_MAX];
	fgets_checked(line,sizeof(line),p_file);

	int n_err = 0;
	read_skip_i(line,&n_err);
	fclose(p_file);

	return n_err;
}

static struct Conv_Order_Data* constructor_Conv_Order_Data
	(const struct Integration_Test_Info*const int_test_info, const int error_type,
	 const struct Simulation*const sim)
{
	const int*const p_range  = get_conv_order_range(int_test_info,'p'),
	         *const ml_range = get_conv_order_range(int_test_info,'m');

	const char*const input_name = compute_error_input_name_root(int_test_info,error_type,false,sim);

	const int n_err = compute_n_err(input_name);
	ptrdiff_t extents[]    = { ml_range[1]+1, p_range[1]+1, n_err, };
	const ptrdiff_t base_e = extents[0]*extents[1];

	char**const var_names = malloc((unsigned)n_err * sizeof *var_names); // moved
	for (int i = 0; i < extents[2]; ++i)
		var_names[i] = malloc(STRLEN_MIN * sizeof *var_names); // moved
	struct Multiarray_i* cases_run   = constructor_zero_Multiarray_i('C',2,&extents[0]); // moved
	struct Multiarray_d* h           = constructor_zero_Multiarray_d('C',2,&extents[0]); // moved
	struct Multiarray_i* ex_ord      = constructor_zero_Multiarray_i('C',3,extents);     // moved
	struct Multiarray_d* l2_err      = constructor_zero_Multiarray_d('C',3,extents);     // moved
	struct Multiarray_d* conv_orders = constructor_zero_Multiarray_d('C',3,extents);     // moved

	bool read_var_names = true;
	for (int p = p_range[0]; p <= p_range[1]; ++p) {
	for (int ml = ml_range[0]; ml <= ml_range[1]; ++ml) {
		const char*const input_name_curr = set_file_name_curr(ADAPT_0,p,ml,true,input_name);

		char line[STRLEN_MAX];
		FILE* p_file = fopen_sp_input_file_unchecked('p',input_name_curr,"txt",0); // closed
		if (!p_file)
			continue;

		skip_lines(p_file,1);
		fgets_checked(line,sizeof(line),p_file);
		int data_i[n_err];
		read_skip_i_1(line,1,data_i,n_err);

		if (read_var_names) {
			read_var_names = false;
			skip_lines(p_file,1);
			fgets_checked(line,sizeof(line),p_file);
			read_skip_c_2(line,1,var_names,n_err);
		} else {
			skip_lines(p_file,2);
		}

		fgets_checked(line,sizeof(line),p_file);

		double data_d[n_err+1];
		read_skip_d_1(line,0,data_d,n_err+1);
		fclose(p_file);

		const ptrdiff_t ind_h = compute_index_sub_container(h->order,0,h->extents,(ptrdiff_t[]){ml,p});
		h->data[ind_h] = data_d[0];
		cases_run->data[ind_h] = 1;

		for (int n = 0; n < n_err; ++n) {
			const ptrdiff_t ind_e = n*base_e + p*extents[0] + ml;
			l2_err->data[ind_e] = data_d[n+1];
			ex_ord->data[ind_e] = data_i[n];
		}
	}}

	for (int p = p_range[0]; p <= p_range[1]; ++p) {
	for (int ml = ml_range[0]+1; ml <= ml_range[1]; ++ml) {
	for (int n = 0; n < n_err; ++n) {
		const ptrdiff_t ind_c_h = p*extents[0] + ml-1,
		                ind_f_h = p*extents[0] + ml,
		                ind_c_e = n*base_e + ind_c_h,
		                ind_f_e = n*base_e + ind_f_h;

		double* data_l2 = l2_err->data,
		      * data_h  = h->data;

		if (data_h[ind_f_h] == 0.0 || data_h[ind_c_h] == 0.0)
			continue;

		conv_orders->data[ind_f_e] =
			log10(data_l2[ind_f_e]/data_l2[ind_c_e])/log10(data_h[ind_f_h]/data_h[ind_c_h]);
	}}}

	struct Conv_Order_Data*const cod = malloc(sizeof *cod);      // free
	cod->var_names   = (const char*const*) var_names;            // free (2 levels)
	cod->cases_run   = (struct const_Multiarray_i*) cases_run;   // destructed
	cod->h           = (struct const_Multiarray_d*) h;           // destructed
	cod->ex_ord      = (struct const_Multiarray_i*) ex_ord;      // destructed
	cod->l2_err      = (struct const_Multiarray_d*) l2_err;      // destructed
	cod->conv_orders = (struct const_Multiarray_d*) conv_orders; // destructed

	return cod;
}

static void destructor_Conv_Order_Data (const struct Conv_Order_Data*const cod)
{
	const ptrdiff_t n_var = cod->l2_err->extents[2];
	for (int i = 0; i < n_var; ++i)
		free((void*)cod->var_names[i]);
	free((void*)cod->var_names);
	destructor_const_Multiarray_i(cod->cases_run);
	destructor_const_Multiarray_d(cod->h);
	destructor_const_Multiarray_i(cod->ex_ord);
	destructor_const_Multiarray_d(cod->l2_err);
	destructor_const_Multiarray_d(cod->conv_orders);
	free((void*)cod);
}

static bool attained_expected_conv_orders
	(const double discount, const struct const_Multiarray_d*const conv_orders,
	 const struct const_Multiarray_i*const exp_orders, const struct Integration_Test_Info* int_test_info)
{
	bool pass = true;

	const ptrdiff_t*const extents_e = conv_orders->extents;
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;
	const int n_err   = (int)extents_e[2];
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

static void output_combined_results
	(const struct Conv_Order_Data*const cod, const struct Integration_Test_Info*const int_test_info,
	 const int error_type, const struct Simulation*const sim)
{
	char error_ext[STRLEN_MIN];
	switch (error_type) {
		case ERROR_STANDARD: strcpy(error_ext,"standard"); break;
		case ERROR_FUNCTIONAL: strcpy(error_ext,"functional"); break;
		default: EXIT_ERROR("Unsupported: %d.\n",error_type); break;
	}

	const char*const output_path = extract_path(compute_error_input_name_root(int_test_info,error_type,false,sim));
	char output_name[STRLEN_MAX];
	sprintf(output_name,"%s%c%s%c%s%s",output_path,'/',"l2errs+convergence",'_',error_ext,".txt");

	const ptrdiff_t n_vars = cod->conv_orders->extents[2];

	FILE* file = fopen_create_dir(output_name); // closed

	fprintf(file,"n_vars %2td\n",n_vars);
	fprintf(file,"ml_max %2d\n",ML_MAX);
	fprintf(file,"p_max  %2d\n",P_MAX);

	struct Matrix_i tmp_i;
	struct Matrix_d tmp_d;

	fprintf(file,"\nCases Run\n");
	set_const_Matrix_from_Multiarray_i((struct const_Matrix_i*)&tmp_i,cod->cases_run,0);
	print_to_file_Matrix_i(file,&tmp_i);

	fprintf(file,"\nMesh Size\n");
	set_const_Matrix_from_Multiarray_d((struct const_Matrix_d*)&tmp_d,cod->h,0);
	print_to_file_Matrix_d(file,&tmp_d);

	fprintf(file,"\nL2 Errors\n\n");
	for (int i = 0; i < n_vars; ++i) {
		fprintf(file,"%s\n",cod->var_names[i]);
		set_const_Matrix_from_Multiarray_d((struct const_Matrix_d*)&tmp_d,cod->l2_err,(ptrdiff_t[]){i});
		print_to_file_Matrix_d(file,&tmp_d);
	}

	fprintf(file,"\nConvergence Orders\n\n");
	for (int i = 0; i < n_vars; ++i) {
		fprintf(file,"%s\n",cod->var_names[i]);
		set_const_Matrix_from_Multiarray_d((struct const_Matrix_d*)&tmp_d,cod->conv_orders,(ptrdiff_t[]){i});
		print_to_file_Matrix_d(file,&tmp_d);
	}

	fclose(file);
}

// Level 2 ********************************************************************************************************** //

static const int* get_conv_order_range (const struct Integration_Test_Info*const int_test_info, const char mp_type)
{
	assert(mp_type == 'm' || mp_type == 'p');

	const char*const conv_study_ext = int_test_info->conv_study_extension;
	if (!conv_study_ext) {
		if      (mp_type == 'p') return int_test_info->p_ref;
		else if (mp_type == 'm') return int_test_info->ml;
	} else {
		static const int max_p_range[]  = { 0, P_MAX, },
		                 max_ml_range[] = { 0, ML_MAX, };
		if      (mp_type == 'p') return max_p_range;
		else if (mp_type == 'm') return max_ml_range;
	}
	EXIT_ERROR("Should not have made it here (mp_type = %c).\n",mp_type);
}

static const char* compute_error_input_name_root
	(const struct Integration_Test_Info*const int_test_info, const int error_type, const bool use_default,
	 const struct Simulation*const sim)
{
	const char*const conv_study_ext = int_test_info->conv_study_extension;
	if (use_default || !conv_study_ext)
		return compute_error_file_name(error_type,sim);
	else {
		static const char* name_part = "../output/results/";

		char l2_spec[STRLEN_MIN];
		switch (error_type) {
			case ERROR_STANDARD: strcpy(l2_spec,""); break;
			case ERROR_FUNCTIONAL: strcpy(l2_spec,"functional_"); break;
			default: EXIT_ERROR("Unsupported: %d\n",error_type); break;
		}

		static char output_name[STRLEN_MAX];
		sprintf(output_name,"%s%s%c%s%s",name_part,conv_study_ext,'/',l2_spec,"l2_errors");
		correct_file_name_ml_p(sim->ml_p_curr[0],sim->ml_p_curr[1],output_name);
		return output_name;
	}
	return NULL;
}
