// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "test_unit_cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "test_base.h"
#include "test_support_matrix.h"
#include "test_support_cubature.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_cubature.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "vector.h"

#include "cubature.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the tensor-product basis functions.
static void test_unit_cubature_tensor_product
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_cubature (struct Test_Info*const test_info)
{
	test_unit_cubature_tensor_product(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for cubature data to be tested for comparison with expected values.
struct Cubature_Data {
	const struct const_Cubature* d1_p3_EQ,  ///< Self-explanatory.
	                           * d2_p4_EQ,  ///< Self-explanatory.
	                           * d3_p3_EQ,  ///< Self-explanatory.
	                           * d1_p3_GLL, ///< Self-explanatory.
	                           * d2_p4_GLL, ///< Self-explanatory.
	                           * d3_p3_GLL, ///< Self-explanatory.
	                           * d1_p3_GL,  ///< Self-explanatory.
	                           * d2_p4_GL,  ///< Self-explanatory.
	                           * d3_p3_GL;  ///< Self-explanatory.
};

/** \brief Constructor for the \ref Cubature_Data.
 *  \return Standard. */
static struct Cubature_Data* constructor_Cubature_Data
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const cub_type  ///< The type of cubature to test.
	);

/// \brief Destructor for the \ref Cubature_Data.
static void destructor_Cubature_Data
	(struct Cubature_Data* cub_data ///< Standard.
	);

static void test_unit_cubature_tensor_product (struct Test_Info*const test_info)
{
	char* test_name = "cubature_tp";
	char test_name_full[STRLEN_MAX];
	strcpy(test_name_full,"Cubature - ");
	strcat(test_name_full,test_name);

	bool pass = true;

	struct Cubature_Data* cub_data_r = constructor_Cubature_Data('r',test_name), // destructed
	                    * cub_data_c = constructor_Cubature_Data('c',test_name); // destructed

	if ((diff_const_Cubature(cub_data_r->d1_p3_EQ, cub_data_c->d1_p3_EQ, EPS)) ||
	    (diff_const_Cubature(cub_data_r->d2_p4_EQ, cub_data_c->d2_p4_EQ, EPS)) ||
	    (diff_const_Cubature(cub_data_r->d3_p3_EQ, cub_data_c->d3_p3_EQ, EPS)) ||
	    (diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,EPS)) ||
	    (diff_const_Cubature(cub_data_r->d2_p4_GLL,cub_data_c->d2_p4_GLL,EPS)) ||
	    (diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,EPS*10)) ||
	    (diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,EPS)))
	{
		test_print_failure(test_info,test_name_full);

		pass = false;

		if (diff_const_Cubature(cub_data_r->d1_p3_EQ,cub_data_c->d1_p3_EQ,EPS))
			print_diff_const_Cubature(cub_data_r->d1_p3_EQ,cub_data_c->d1_p3_EQ,EPS);
		if (diff_const_Cubature(cub_data_r->d2_p4_EQ,cub_data_c->d2_p4_EQ,EPS))
			print_diff_const_Cubature(cub_data_r->d2_p4_EQ,cub_data_c->d2_p4_EQ,EPS);
		if (diff_const_Cubature(cub_data_r->d3_p3_EQ,cub_data_c->d3_p3_EQ,EPS))
			print_diff_const_Cubature(cub_data_r->d3_p3_EQ,cub_data_c->d3_p3_EQ,EPS);

		if (diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,EPS))
			print_diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,EPS);
		if (diff_const_Cubature(cub_data_r->d2_p4_GLL,cub_data_c->d2_p4_GLL,EPS))
			print_diff_const_Cubature(cub_data_r->d2_p4_GLL,cub_data_c->d2_p4_GLL,EPS);
		if (diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,EPS))
			print_diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,EPS);

		if (diff_const_Cubature(cub_data_r->d1_p3_GL,cub_data_c->d1_p3_GL,EPS))
			print_diff_const_Cubature(cub_data_r->d1_p3_GL,cub_data_c->d1_p3_GL,EPS);
	}

	destructor_Cubature_Data(cub_data_r);
	destructor_Cubature_Data(cub_data_c);

	test_increment_and_print(test_info,pass,test_name_full);
}

// Level 1 ********************************************************************************************************** //

static struct Cubature_Data* constructor_Cubature_Data (const char eval_type, const char*const cub_type)
{
	struct Cubature_Data* cub_data = calloc(1,sizeof *cub_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"cubature/");
		strcat(file_name_part,cub_type);
		const char*const file_name_full = constructor_file_name_unit(file_name_part); // free

		cub_data->d1_p3_EQ  = constructor_file_name_const_Cubature("Cubature_d1_p3_EQ",file_name_full);   // keep
		cub_data->d2_p4_EQ  = constructor_file_name_const_Cubature("Cubature_d2_p4_EQ",file_name_full);   // keep
		cub_data->d3_p3_EQ  = constructor_file_name_const_Cubature("Cubature_d3_p3_EQ",file_name_full);   // keep
		cub_data->d1_p3_GLL = constructor_file_name_const_Cubature("Cubature_d1_p3_GLoL",file_name_full); // keep
		cub_data->d2_p4_GLL = constructor_file_name_const_Cubature("Cubature_d2_p4_GLoL",file_name_full); // keep
		cub_data->d3_p3_GLL = constructor_file_name_const_Cubature("Cubature_d3_p3_GLoL",file_name_full); // keep
		cub_data->d1_p3_GL  = constructor_file_name_const_Cubature("Cubature_d1_p3_GLe",file_name_full);  // keep

		free((void*)file_name_full);
	} else if (eval_type == 'c') {
		cub_data->d1_p3_EQ  = constructor_const_Cubature_TP(1,3,CUB_EQ);  // keep
		cub_data->d2_p4_EQ  = constructor_const_Cubature_TP(2,4,CUB_EQ);  // keep
		cub_data->d3_p3_EQ  = constructor_const_Cubature_TP(3,3,CUB_EQ);  // keep
		cub_data->d1_p3_GLL = constructor_const_Cubature_TP(1,3,CUB_GLL); // keep
		cub_data->d2_p4_GLL = constructor_const_Cubature_TP(2,4,CUB_GLL); // keep
		cub_data->d3_p3_GLL = constructor_const_Cubature_TP(3,3,CUB_GLL); // keep
		cub_data->d1_p3_GL  = constructor_const_Cubature_TP(1,3,CUB_GL);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return cub_data;
}

static void destructor_Cubature_Data (struct Cubature_Data* cub_data)
{
	destructor_const_Cubature(cub_data->d1_p3_EQ);
	destructor_const_Cubature(cub_data->d2_p4_EQ);
	destructor_const_Cubature(cub_data->d3_p3_EQ);
	destructor_const_Cubature(cub_data->d1_p3_GLL);
	destructor_const_Cubature(cub_data->d2_p4_GLL);
	destructor_const_Cubature(cub_data->d3_p3_GLL);
	destructor_const_Cubature(cub_data->d1_p3_GL);
	free(cub_data);
}
