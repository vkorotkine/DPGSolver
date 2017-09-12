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
#include "test_support_cubature.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_cubature.h"
#include "definitions_tol.h"

#include "cubature.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the tensor-product cubature.
static void test_unit_cubature_tensor_product
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the simplex cubature.
static void test_unit_cubature_simplex
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the pyramid cubature.
static void test_unit_cubature_pyramid
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

void test_unit_cubature (struct Test_Info*const test_info)
{
	test_unit_cubature_tensor_product(test_info);
	test_unit_cubature_simplex(test_info);
	test_unit_cubature_pyramid(test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for tensor-product cubature data to be tested for comparison with expected values.
struct Cubature_Data_TP {
	const struct const_Cubature* d1_p3_EQ,  ///< See \ref definitions_cubature.h.
	                           * d2_p4_EQ,  ///< See \ref definitions_cubature.h.
	                           * d3_p3_EQ,  ///< See \ref definitions_cubature.h.
	                           * d1_p3_GLL, ///< See \ref definitions_cubature.h.
	                           * d2_p4_GLL, ///< See \ref definitions_cubature.h.
	                           * d3_p3_GLL, ///< See \ref definitions_cubature.h.
	                           * d1_p3_GL,  ///< See \ref definitions_cubature.h.
	                           * d2_p4_GL,  ///< See \ref definitions_cubature.h.
	                           * d3_p3_GL;  ///< See \ref definitions_cubature.h.
};

/// Container for simplex cubature data to be tested for comparison with expected values.
struct Cubature_Data_SI {
	const struct const_Cubature* d2_p3_AO,  ///< See \ref definitions_cubature.h.
	                           * d2_p3_WSH, ///< See \ref definitions_cubature.h.
	                           * d2_p3_EQ,  ///< See \ref definitions_cubature.h.
	                           * d2_p6_WV,  ///< See \ref definitions_cubature.h.
	                           * d3_p2_AO,  ///< See \ref definitions_cubature.h.
	                           * d3_p2_WSH, ///< See \ref definitions_cubature.h.
	                           * d3_p4_WV;  ///< See \ref definitions_cubature.h.
};

/// Container for pyramid cubature data to be tested for comparison with expected values.
struct Cubature_Data_PYR {
	const struct const_Cubature* d3_p2_GL,  ///< See \ref definitions_cubature.h.
	                           * d3_p3_GLL, ///< See \ref definitions_cubature.h.
	                           * d3_p4_WV;  ///< See \ref definitions_cubature.h.
};

/** \brief Constructor for the \ref Cubature_Data_TP.
 *  \return Standard. */
static struct Cubature_Data_TP* constructor_Cubature_Data_TP
	(const char eval_type,     ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const cub_type ///< The type of cubature to test.
	);

/// \brief Destructor for the \ref Cubature_Data_TP.
static void destructor_Cubature_Data_TP
	(struct Cubature_Data_TP* cub_data ///< Standard.
	);

/** \brief Constructor for the \ref Cubature_Data_SI.
 *  \return Standard. */
static struct Cubature_Data_SI* constructor_Cubature_Data_SI
	(const char eval_type,     ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const cub_type ///< The type of cubature to test.
	);

/// \brief Destructor for the \ref Cubature_Data_SI.
static void destructor_Cubature_Data_SI
	(struct Cubature_Data_SI* cub_data ///< Standard.
	);

/** \brief Constructor for the \ref Cubature_Data_PYR.
 *  \return Standard. */
static struct Cubature_Data_PYR* constructor_Cubature_Data_PYR
	(const char eval_type,     ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const cub_type ///< The type of cubature to test.
	);

/// \brief Destructor for the \ref Cubature_Data_PYR.
static void destructor_Cubature_Data_PYR
	(struct Cubature_Data_PYR* cub_data ///< Standard.
	);

static void test_unit_cubature_tensor_product (struct Test_Info*const test_info)
{
	char* test_name = "cubature_tp";
	char test_name_full[STRLEN_MAX];
	strcpy(test_name_full,"Cubature - ");
	strcat(test_name_full,test_name);

	bool pass = true;

	struct Cubature_Data_TP* cub_data_r = constructor_Cubature_Data_TP('r',test_name), // destructed
	                       * cub_data_c = constructor_Cubature_Data_TP('c',test_name); // destructed

	double tol[]       = { EPS, EPS, EPS, EPS, EPS, 10*EPS, EPS, 10*EPS, 10*EPS };
	bool differences[] =
		{ diff_const_Cubature(cub_data_r->d1_p3_EQ, cub_data_c->d1_p3_EQ, tol[0]),
		  diff_const_Cubature(cub_data_r->d2_p4_EQ, cub_data_c->d2_p4_EQ, tol[1]),
		  diff_const_Cubature(cub_data_r->d3_p3_EQ, cub_data_c->d3_p3_EQ, tol[2]),
		  diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,tol[3]),
		  diff_const_Cubature(cub_data_r->d2_p4_GLL,cub_data_c->d2_p4_GLL,tol[4]),
		  diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,tol[5]),
		  diff_const_Cubature(cub_data_r->d1_p3_GL, cub_data_c->d1_p3_GL, tol[6]),
		  diff_const_Cubature(cub_data_r->d2_p4_GL, cub_data_c->d2_p4_GL, tol[7]),
		  diff_const_Cubature(cub_data_r->d3_p3_GL, cub_data_c->d3_p3_GL, tol[8]),
		};

	const int n_diff = sizeof(differences)/sizeof(*differences);

	bool diff = false;
	for (int i = 0; i < n_diff; ++i) {
		if (differences[i]) {
			diff = true;
			break;
		}
	}

	if (diff) {
		test_print_failure(test_info,test_name_full);

		pass = false;

		if (differences[0])
			print_diff_const_Cubature(cub_data_r->d1_p3_EQ,cub_data_c->d1_p3_EQ,tol[0]);
		if (differences[1])
			print_diff_const_Cubature(cub_data_r->d2_p4_EQ,cub_data_c->d2_p4_EQ,tol[1]);
		if (differences[2])
			print_diff_const_Cubature(cub_data_r->d3_p3_EQ,cub_data_c->d3_p3_EQ,tol[2]);

		if (differences[3])
			print_diff_const_Cubature(cub_data_r->d1_p3_GLL,cub_data_c->d1_p3_GLL,tol[3]);
		if (differences[4])
			print_diff_const_Cubature(cub_data_r->d2_p4_GLL,cub_data_c->d2_p4_GLL,tol[4]);
		if (differences[5])
			print_diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,tol[5]);

		if (differences[6])
			print_diff_const_Cubature(cub_data_r->d1_p3_GL,cub_data_c->d1_p3_GL,tol[6]);
		if (differences[7])
			print_diff_const_Cubature(cub_data_r->d2_p4_GL,cub_data_c->d2_p4_GL,tol[7]);
		if (differences[8])
			print_diff_const_Cubature(cub_data_r->d3_p3_GL,cub_data_c->d3_p3_GL,tol[8]);
	}

	destructor_Cubature_Data_TP(cub_data_r);
	destructor_Cubature_Data_TP(cub_data_c);

	test_increment_and_print(test_info,pass,test_name_full);
}

static void test_unit_cubature_simplex (struct Test_Info*const test_info)
{
	char* test_name = "cubature_si";
	char test_name_full[STRLEN_MAX];
	strcpy(test_name_full,"Cubature - ");
	strcat(test_name_full,test_name);

	bool pass = true;

	struct Cubature_Data_SI* cub_data_r = constructor_Cubature_Data_SI('r',test_name), // destructed
	                       * cub_data_c = constructor_Cubature_Data_SI('c',test_name); // destructed

	double tol[]       = { 10*EPS, 10*EPS, 10*EPS, 1e2*EPS, EPS, 10*EPS, 10*EPS, };
	bool differences[] =
		{ diff_const_Cubature(cub_data_r->d2_p3_EQ, cub_data_c->d2_p3_EQ, tol[0]),
		  diff_const_Cubature(cub_data_r->d2_p3_AO, cub_data_c->d2_p3_AO, tol[1]),
		  diff_const_Cubature(cub_data_r->d2_p3_WSH,cub_data_c->d2_p3_WSH,tol[2]),
		  diff_const_Cubature(cub_data_r->d2_p6_WV, cub_data_c->d2_p6_WV, tol[3]),
		  diff_const_Cubature(cub_data_r->d3_p2_AO, cub_data_c->d3_p2_AO, tol[4]),
		  diff_const_Cubature(cub_data_r->d3_p2_WSH,cub_data_c->d3_p2_WSH,tol[5]),
		  diff_const_Cubature(cub_data_r->d3_p4_WV, cub_data_c->d3_p4_WV, tol[6]),
		};

	const int n_diff = sizeof(differences)/sizeof(*differences);

	bool diff = false;
	for (int i = 0; i < n_diff; ++i) {
		if (differences[i]) {
			diff = true;
			break;
		}
	}

	if (diff) {
		test_print_failure(test_info,test_name_full);

		pass = false;

		if (differences[0]) print_diff_const_Cubature(cub_data_r->d2_p3_EQ, cub_data_c->d2_p3_EQ, tol[0]);
		if (differences[1]) print_diff_const_Cubature(cub_data_r->d2_p3_AO, cub_data_c->d2_p3_AO, tol[1]);
		if (differences[2]) print_diff_const_Cubature(cub_data_r->d2_p3_WSH,cub_data_c->d2_p3_WSH,tol[2]);
		if (differences[3]) print_diff_const_Cubature(cub_data_r->d2_p6_WV, cub_data_c->d2_p6_WV, tol[3]);

		if (differences[4]) print_diff_const_Cubature(cub_data_r->d3_p2_AO, cub_data_c->d3_p2_AO, tol[4]);
		if (differences[5]) print_diff_const_Cubature(cub_data_r->d3_p2_WSH,cub_data_c->d3_p2_WSH,tol[5]);
		if (differences[6]) print_diff_const_Cubature(cub_data_r->d3_p4_WV, cub_data_c->d3_p4_WV, tol[6]);
	}

	destructor_Cubature_Data_SI(cub_data_r);
	destructor_Cubature_Data_SI(cub_data_c);

	test_increment_and_print(test_info,pass,test_name_full);
}

static void test_unit_cubature_pyramid (struct Test_Info*const test_info)
{
	char* test_name = "cubature_pyr";
	char test_name_full[STRLEN_MAX];
	strcpy(test_name_full,"Cubature - ");
	strcat(test_name_full,test_name);

	bool pass = true;

	struct Cubature_Data_PYR* cub_data_r = constructor_Cubature_Data_PYR('r',test_name), // destructed
	                        * cub_data_c = constructor_Cubature_Data_PYR('c',test_name); // destructed

	double tol[]       = { 10*EPS, 10*EPS, 10*EPS, };
	bool differences[] =
		{ diff_const_Cubature(cub_data_r->d3_p2_GL, cub_data_c->d3_p2_GL, tol[0]),
		  diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,tol[1]),
		  diff_const_Cubature(cub_data_r->d3_p4_WV, cub_data_c->d3_p4_WV, tol[2]),
		};

	test_print_warning(test_info,"Not all PYR nodes being tested.");

	const int n_diff = sizeof(differences)/sizeof(*differences);

	bool diff = false;
	for (int i = 0; i < n_diff; ++i) {
		if (differences[i]) {
			diff = true;
			break;
		}
	}

	if (diff) {
		test_print_failure(test_info,test_name_full);

		pass = false;

		if (differences[0]) print_diff_const_Cubature(cub_data_r->d3_p2_GL, cub_data_c->d3_p2_GL, tol[0]);
		if (differences[1]) print_diff_const_Cubature(cub_data_r->d3_p3_GLL,cub_data_c->d3_p3_GLL,tol[1]);
		if (differences[2]) print_diff_const_Cubature(cub_data_r->d3_p4_WV, cub_data_c->d3_p4_WV, tol[2]);
	}

	destructor_Cubature_Data_PYR(cub_data_r);
	destructor_Cubature_Data_PYR(cub_data_c);

	test_increment_and_print(test_info,pass,test_name_full);
}

// Level 1 ********************************************************************************************************** //

static struct Cubature_Data_TP* constructor_Cubature_Data_TP (const char eval_type, const char*const cub_type)
{
	struct Cubature_Data_TP* cub_data = calloc(1,sizeof *cub_data); // returned
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
		cub_data->d2_p4_GL  = constructor_file_name_const_Cubature("Cubature_d2_p4_GLe",file_name_full);  // keep
		cub_data->d3_p3_GL  = constructor_file_name_const_Cubature("Cubature_d3_p3_GLe",file_name_full);  // keep

		free((void*)file_name_full);
	} else if (eval_type == 'c') {
		cub_data->d1_p3_EQ  = constructor_const_Cubature_tp(1,3,CUB_EQ);  // keep
		cub_data->d2_p4_EQ  = constructor_const_Cubature_tp(2,4,CUB_EQ);  // keep
		cub_data->d3_p3_EQ  = constructor_const_Cubature_tp(3,3,CUB_EQ);  // keep
		cub_data->d1_p3_GLL = constructor_const_Cubature_tp(1,3,CUB_GLL); // keep
		cub_data->d2_p4_GLL = constructor_const_Cubature_tp(2,4,CUB_GLL); // keep
		cub_data->d3_p3_GLL = constructor_const_Cubature_tp(3,3,CUB_GLL); // keep
		cub_data->d1_p3_GL  = constructor_const_Cubature_tp(1,3,CUB_GL);  // keep
		cub_data->d2_p4_GL  = constructor_const_Cubature_tp(2,4,CUB_GL);  // keep
		cub_data->d3_p3_GL  = constructor_const_Cubature_tp(3,3,CUB_GL);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return cub_data;
}

static void destructor_Cubature_Data_TP (struct Cubature_Data_TP* cub_data)
{
	destructor_const_Cubature(cub_data->d1_p3_EQ);
	destructor_const_Cubature(cub_data->d2_p4_EQ);
	destructor_const_Cubature(cub_data->d3_p3_EQ);
	destructor_const_Cubature(cub_data->d1_p3_GLL);
	destructor_const_Cubature(cub_data->d2_p4_GLL);
	destructor_const_Cubature(cub_data->d3_p3_GLL);
	destructor_const_Cubature(cub_data->d1_p3_GL);
	destructor_const_Cubature(cub_data->d2_p4_GL);
	destructor_const_Cubature(cub_data->d3_p3_GL);

	free(cub_data);
}

static struct Cubature_Data_SI* constructor_Cubature_Data_SI (const char eval_type, const char*const cub_type)
{
	struct Cubature_Data_SI* cub_data = calloc(1,sizeof *cub_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"cubature/");
		strcat(file_name_part,cub_type);
		const char*const file_name_full = constructor_file_name_unit(file_name_part); // free

		cub_data->d2_p3_EQ  = constructor_file_name_const_Cubature("Cubature_d2_p3_EQ",file_name_full);  // keep
		cub_data->d2_p3_AO  = constructor_file_name_const_Cubature("Cubature_d2_p3_AO",file_name_full);  // keep
		cub_data->d2_p3_WSH = constructor_file_name_const_Cubature("Cubature_d2_p3_WSH",file_name_full); // keep
		cub_data->d2_p6_WV  = constructor_file_name_const_Cubature("Cubature_d2_p6_WV",file_name_full);  // keep
		cub_data->d3_p2_AO  = constructor_file_name_const_Cubature("Cubature_d3_p2_AO",file_name_full);  // keep
		cub_data->d3_p2_WSH = constructor_file_name_const_Cubature("Cubature_d3_p2_WSH",file_name_full); // keep
		cub_data->d3_p4_WV  = constructor_file_name_const_Cubature("Cubature_d3_p4_WV",file_name_full);  // keep

		free((void*)file_name_full);
	} else if (eval_type == 'c') {
		cub_data->d2_p3_EQ  = constructor_const_Cubature_si(2,3,CUB_EQ);  // keep
		cub_data->d2_p3_AO  = constructor_const_Cubature_si(2,3,CUB_AO);  // keep
		cub_data->d2_p3_WSH = constructor_const_Cubature_si(2,3,CUB_WSH); // keep
		cub_data->d2_p6_WV  = constructor_const_Cubature_si(2,6,CUB_WV);  // keep
		cub_data->d3_p2_AO  = constructor_const_Cubature_si(3,2,CUB_AO);  // keep
		cub_data->d3_p2_WSH = constructor_const_Cubature_si(3,2,CUB_WSH); // keep
		cub_data->d3_p4_WV  = constructor_const_Cubature_si(3,4,CUB_WV);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return cub_data;
}

static void destructor_Cubature_Data_SI (struct Cubature_Data_SI* cub_data)
{
	destructor_const_Cubature(cub_data->d2_p3_EQ);
	destructor_const_Cubature(cub_data->d2_p3_AO);
	destructor_const_Cubature(cub_data->d2_p3_WSH);
	destructor_const_Cubature(cub_data->d2_p6_WV);
	destructor_const_Cubature(cub_data->d3_p2_AO);
	destructor_const_Cubature(cub_data->d3_p2_WSH);
	destructor_const_Cubature(cub_data->d3_p4_WV);

	free(cub_data);
}

static struct Cubature_Data_PYR* constructor_Cubature_Data_PYR (const char eval_type, const char*const cub_type)
{
	struct Cubature_Data_PYR* cub_data = calloc(1,sizeof *cub_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"cubature/");
		strcat(file_name_part,cub_type);
		const char*const file_name_full = constructor_file_name_unit(file_name_part); // free

		cub_data->d3_p2_GL  = constructor_file_name_const_Cubature("Cubature_d3_p2_GLe", file_name_full); // keep
		cub_data->d3_p3_GLL = constructor_file_name_const_Cubature("Cubature_d3_p3_GLoL",file_name_full); // keep
		cub_data->d3_p4_WV  = constructor_file_name_const_Cubature("Cubature_d3_p4_WV",  file_name_full); // keep

		free((void*)file_name_full);
	} else if (eval_type == 'c') {
		cub_data->d3_p2_GL  = constructor_const_Cubature_pyr(3,2,CUB_GL);  // keep
		cub_data->d3_p3_GLL = constructor_const_Cubature_pyr(3,3,CUB_GLL); // keep
		cub_data->d3_p4_WV  = constructor_const_Cubature_pyr(3,4,CUB_WV);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return cub_data;
}

static void destructor_Cubature_Data_PYR (struct Cubature_Data_PYR* cub_data)
{
	destructor_const_Cubature(cub_data->d3_p2_GL);
	destructor_const_Cubature(cub_data->d3_p3_GLL);
	destructor_const_Cubature(cub_data->d3_p4_WV);

	free(cub_data);
}
