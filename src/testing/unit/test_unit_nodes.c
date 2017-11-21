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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "test_base.h"
#include "test_support.h"
#include "test_support_nodes.h"
#include "test_support_multiarray.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_nodes.h"
#include "definitions_tol.h"
#include "definitions_elements.h"

#include "multiarray.h"

#include "nodes.h"
#include "nodes_plotting.h"
#include "nodes_correspondence.h"

// Static function declarations ************************************************************************************* //

/// \brief Provides unit tests for the tensor-product nodes.
static void test_unit_nodes_tensor_product
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the simplex nodes.
static void test_unit_nodes_simplex
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the pyramid nodes.
static void test_unit_nodes_pyramid
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the plotting nodes for all element types.
static void test_unit_nodes_plotting
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/// \brief Provides unit tests for the node correspondence for all permutations of adjacent faces.
static void test_unit_nodes_face_correspondence
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the nodes (\ref test_unit_nodes.c).
 *  \return 0 on success. */
int main
	(int nargc,  ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(nargc == 2,"Invalid number of input arguments");
	const char* test_name = argv[1];

	struct Test_Info test_info = { .n_warn = 0, };
	if (strcmp(test_name,"tp") == 0)
		test_unit_nodes_tensor_product(&test_info);
	else if (strcmp(test_name,"si") == 0)
		test_unit_nodes_simplex(&test_info);
	else if (strcmp(test_name,"pyr") == 0)
		test_unit_nodes_pyramid(&test_info);
	else if (strcmp(test_name,"plotting") == 0)
		test_unit_nodes_plotting(&test_info);
	else if (strcmp(test_name,"face_correspondence") == 0)
		test_unit_nodes_face_correspondence(&test_info);
	else
		EXIT_ERROR("Invalid test name: %s\n",test_name);
	output_warning_count(&test_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// Container for tensor-product node data to be tested for comparison with expected values.
struct Nodes_Data_TP {
	const struct const_Nodes* d1_p3_EQ,  ///< See \ref definitions_nodes.h.
	                        * d2_p4_EQ,  ///< See \ref definitions_nodes.h.
	                        * d3_p3_EQ,  ///< See \ref definitions_nodes.h.
	                        * d1_p3_GLL, ///< See \ref definitions_nodes.h.
	                        * d2_p4_GLL, ///< See \ref definitions_nodes.h.
	                        * d3_p3_GLL, ///< See \ref definitions_nodes.h.
	                        * d1_p3_GL,  ///< See \ref definitions_nodes.h.
	                        * d2_p4_GL,  ///< See \ref definitions_nodes.h.
	                        * d3_p3_GL;  ///< See \ref definitions_nodes.h.
};

/// Container for simplex node data to be tested for comparison with expected values.
struct Nodes_Data_SI {
	const struct const_Nodes* d2_p3_AO,  ///< See \ref definitions_nodes.h.
	                        * d2_p3_WSH, ///< See \ref definitions_nodes.h.
	                        * d2_p3_EQ,  ///< See \ref definitions_nodes.h.
	                        * d2_p6_WV,  ///< See \ref definitions_nodes.h.
	                        * d3_p2_AO,  ///< See \ref definitions_nodes.h.
	                        * d3_p2_WSH, ///< See \ref definitions_nodes.h.
	                        * d3_p4_WV;  ///< See \ref definitions_nodes.h.
};

/// Container for pyramid node data to be tested for comparison with expected values.
struct Nodes_Data_PYR {
	const struct const_Nodes* d3_p2_GL,  ///< See \ref definitions_nodes.h.
	                        * d3_p3_GLL, ///< See \ref definitions_nodes.h.
	                        * d3_p4_WV;  ///< See \ref definitions_nodes.h.
};

/// Container for plotting node data to be tested for comparison with expected values.
struct Plotting_Nodes_Data {
	const struct const_Plotting_Nodes* plt_line,  ///< Plotting nodes for the LINE  element.
	                                 * plt_tri,   ///< Plotting nodes for the TRI   element.
	                                 * plt_quad,  ///< Plotting nodes for the QUAD  element.
	                                 * plt_tet,   ///< Plotting nodes for the TET   element.
	                                 * plt_hex,   ///< Plotting nodes for the HEX   element.
	                                 * plt_wedge, ///< Plotting nodes for the WEDGE element.
	                                 * plt_pyr;   ///< Plotting nodes for the PYR   element.
};

/// Container for face correspondence node data to be tested for comparison with expected values.
struct Nodes_FC_Data {
	const struct const_Multiarray_Vector_i* point_N, ///< Vectors of correspondence indices for a N- node point face.
	                                      * line_3,  ///< Vectors of correspondence indices for a 3- node line  face.
	                                      * line_4,  ///< Vectors of correspondence indices for a 4- node line  face.
	                                      * quad_9,  ///< Vectors of correspondence indices for a 9- node quad  face.
	                                      * quad_16, ///< Vectors of correspondence indices for a 16-node quad  face.
	                                      * tri_6,   ///< Vectors of correspondence indices for a 6- node tri   face.
	                                      * tri_10;  ///< Vectors of correspondence indices for a 10-node tri   face.
};

/** \brief Constructor for the \ref Nodes_Data_TP.
 *  \return Standard. */
static struct Nodes_Data_TP* constructor_Nodes_Data_TP
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const node_type ///< The type of nodes to test.
	);

/// \brief Destructor for the \ref Nodes_Data_TP.
static void destructor_Nodes_Data_TP
	(struct Nodes_Data_TP* nodes_data ///< Standard.
	);

/** \brief Constructor for the \ref Nodes_Data_SI.
 *  \return Standard. */
static struct Nodes_Data_SI* constructor_Nodes_Data_SI
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const node_type ///< The type of nodes to test.
	);

/// \brief Destructor for the \ref Nodes_Data_SI.
static void destructor_Nodes_Data_SI
	(struct Nodes_Data_SI* nodes_data ///< Standard.
	);

/** \brief Constructor for the \ref Nodes_Data_PYR.
 *  \return Standard. */
static struct Nodes_Data_PYR* constructor_Nodes_Data_PYR
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const node_type ///< The type of nodes to test.
	);

/// \brief Destructor for the \ref Nodes_Data_PYR.
static void destructor_Nodes_Data_PYR
	(struct Nodes_Data_PYR* nodes_data ///< Standard.
	);

/** \brief Constructor for the \ref Plotting_Nodes_Data.
 *  \return Standard. */
static struct Plotting_Nodes_Data* constructor_Plotting_Nodes_Data
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const node_type ///< The type of nodes to test.
	);

/// \brief Destructor for the \ref Plotting_Nodes_Data.
static void destructor_Plotting_Nodes_Data
	(struct Plotting_Nodes_Data* p_nodes_data ///< Standard.
	);

/** \brief Constructor for the \ref Nodes_FC_Data container.
 *  \return Standard. */
static struct Nodes_FC_Data* constructor_Nodes_FC_Data
	(const char eval_type,      ///< Method to use to obtain the data. Options: 'r'ead, 'c'ompute
	 const char*const node_type ///< The type of nodes to test.
	);

/// \brief Destructor for the \ref Nodes_FC_Data.
static void destructor_Nodes_FC_Data
	(struct Nodes_FC_Data* nodes_fc_data ///< Standard.
	);

static void test_unit_nodes_tensor_product (struct Test_Info*const test_info)
{
	char* test_name = "nodes_tp";
	bool pass = true;

	struct Nodes_Data_TP* nodes_data_r = constructor_Nodes_Data_TP('r',test_name), // destructed
	                    * nodes_data_c = constructor_Nodes_Data_TP('c',test_name); // destructed

	double tol[]       = { EPS, EPS, EPS, EPS, EPS, 10*EPS, EPS, 10*EPS, 10*EPS };
	bool differences[] =
		{ diff_const_Nodes(nodes_data_r->d1_p3_EQ, nodes_data_c->d1_p3_EQ, tol[0]),
		  diff_const_Nodes(nodes_data_r->d2_p4_EQ, nodes_data_c->d2_p4_EQ, tol[1]),
		  diff_const_Nodes(nodes_data_r->d3_p3_EQ, nodes_data_c->d3_p3_EQ, tol[2]),
		  diff_const_Nodes(nodes_data_r->d1_p3_GLL,nodes_data_c->d1_p3_GLL,tol[3]),
		  diff_const_Nodes(nodes_data_r->d2_p4_GLL,nodes_data_c->d2_p4_GLL,tol[4]),
		  diff_const_Nodes(nodes_data_r->d3_p3_GLL,nodes_data_c->d3_p3_GLL,tol[5]),
		  diff_const_Nodes(nodes_data_r->d1_p3_GL, nodes_data_c->d1_p3_GL, tol[6]),
		  diff_const_Nodes(nodes_data_r->d2_p4_GL, nodes_data_c->d2_p4_GL, tol[7]),
		  diff_const_Nodes(nodes_data_r->d3_p3_GL, nodes_data_c->d3_p3_GL, tol[8]),
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
		pass = false;

		if (differences[0])
			print_diff_const_Nodes(nodes_data_r->d1_p3_EQ,nodes_data_c->d1_p3_EQ,tol[0]);
		if (differences[1])
			print_diff_const_Nodes(nodes_data_r->d2_p4_EQ,nodes_data_c->d2_p4_EQ,tol[1]);
		if (differences[2])
			print_diff_const_Nodes(nodes_data_r->d3_p3_EQ,nodes_data_c->d3_p3_EQ,tol[2]);

		if (differences[3])
			print_diff_const_Nodes(nodes_data_r->d1_p3_GLL,nodes_data_c->d1_p3_GLL,tol[3]);
		if (differences[4])
			print_diff_const_Nodes(nodes_data_r->d2_p4_GLL,nodes_data_c->d2_p4_GLL,tol[4]);
		if (differences[5])
			print_diff_const_Nodes(nodes_data_r->d3_p3_GLL,nodes_data_c->d3_p3_GLL,tol[5]);

		if (differences[6])
			print_diff_const_Nodes(nodes_data_r->d1_p3_GL,nodes_data_c->d1_p3_GL,tol[6]);
		if (differences[7])
			print_diff_const_Nodes(nodes_data_r->d2_p4_GL,nodes_data_c->d2_p4_GL,tol[7]);
		if (differences[8])
			print_diff_const_Nodes(nodes_data_r->d3_p3_GL,nodes_data_c->d3_p3_GL,tol[8]);
	}

	destructor_Nodes_Data_TP(nodes_data_r);
	destructor_Nodes_Data_TP(nodes_data_c);

	test_print_warning(test_info,"Cubature strengths not currently being tested.");
	assert_condition(pass);
}

static void test_unit_nodes_simplex (struct Test_Info*const test_info)
{
	char* test_name = "nodes_si";
	bool pass = true;

	struct Nodes_Data_SI* nodes_data_r = constructor_Nodes_Data_SI('r',test_name), // destructed
	                    * nodes_data_c = constructor_Nodes_Data_SI('c',test_name); // destructed

	double tol[]       = { 10*EPS, 10*EPS, 10*EPS, 1e2*EPS, EPS, 10*EPS, 10*EPS, };
	bool differences[] =
		{ diff_const_Nodes(nodes_data_r->d2_p3_EQ, nodes_data_c->d2_p3_EQ, tol[0]),
		  diff_const_Nodes(nodes_data_r->d2_p3_AO, nodes_data_c->d2_p3_AO, tol[1]),
		  diff_const_Nodes(nodes_data_r->d2_p3_WSH,nodes_data_c->d2_p3_WSH,tol[2]),
		  diff_const_Nodes(nodes_data_r->d2_p6_WV, nodes_data_c->d2_p6_WV, tol[3]),
		  diff_const_Nodes(nodes_data_r->d3_p2_AO, nodes_data_c->d3_p2_AO, tol[4]),
		  diff_const_Nodes(nodes_data_r->d3_p2_WSH,nodes_data_c->d3_p2_WSH,tol[5]),
		  diff_const_Nodes(nodes_data_r->d3_p4_WV, nodes_data_c->d3_p4_WV, tol[6]),
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
		pass = false;

		if (differences[0]) print_diff_const_Nodes(nodes_data_r->d2_p3_EQ, nodes_data_c->d2_p3_EQ, tol[0]);
		if (differences[1]) print_diff_const_Nodes(nodes_data_r->d2_p3_AO, nodes_data_c->d2_p3_AO, tol[1]);
		if (differences[2]) print_diff_const_Nodes(nodes_data_r->d2_p3_WSH,nodes_data_c->d2_p3_WSH,tol[2]);
		if (differences[3]) print_diff_const_Nodes(nodes_data_r->d2_p6_WV, nodes_data_c->d2_p6_WV, tol[3]);

		if (differences[4]) print_diff_const_Nodes(nodes_data_r->d3_p2_AO, nodes_data_c->d3_p2_AO, tol[4]);
		if (differences[5]) print_diff_const_Nodes(nodes_data_r->d3_p2_WSH,nodes_data_c->d3_p2_WSH,tol[5]);
		if (differences[6]) print_diff_const_Nodes(nodes_data_r->d3_p4_WV, nodes_data_c->d3_p4_WV, tol[6]);
	}

	destructor_Nodes_Data_SI(nodes_data_r);
	destructor_Nodes_Data_SI(nodes_data_c);

	test_print_warning(test_info,"Cubature strengths not currently being tested.");
	assert_condition(pass);
}

static void test_unit_nodes_pyramid (struct Test_Info*const test_info)
{
	char* test_name = "nodes_pyr";
	bool pass = true;

	struct Nodes_Data_PYR* nodes_data_r = constructor_Nodes_Data_PYR('r',test_name), // destructed
	                     * nodes_data_c = constructor_Nodes_Data_PYR('c',test_name); // destructed

	double tol[]       = { 10*EPS, 10*EPS, 10*EPS, };
	bool differences[] =
		{ diff_const_Nodes(nodes_data_r->d3_p2_GL, nodes_data_c->d3_p2_GL, tol[0]),
		  diff_const_Nodes(nodes_data_r->d3_p3_GLL,nodes_data_c->d3_p3_GLL,tol[1]),
		  diff_const_Nodes(nodes_data_r->d3_p4_WV, nodes_data_c->d3_p4_WV, tol[2]),
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
		pass = false;

		if (differences[0]) print_diff_const_Nodes(nodes_data_r->d3_p2_GL, nodes_data_c->d3_p2_GL, tol[0]);
		if (differences[1]) print_diff_const_Nodes(nodes_data_r->d3_p3_GLL,nodes_data_c->d3_p3_GLL,tol[1]);
		if (differences[2]) print_diff_const_Nodes(nodes_data_r->d3_p4_WV, nodes_data_c->d3_p4_WV, tol[2]);
	}

	destructor_Nodes_Data_PYR(nodes_data_r);
	destructor_Nodes_Data_PYR(nodes_data_c);

	test_print_warning(test_info,"Cubature strengths not currently being tested.");
	assert_condition(pass);
}

static void test_unit_nodes_plotting (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	char* test_name = "nodes_plotting";
	bool pass = true;

	struct Plotting_Nodes_Data* p_nodes_data_r = constructor_Plotting_Nodes_Data('r',test_name), // destructed
	                          * p_nodes_data_c = constructor_Plotting_Nodes_Data('c',test_name); // destructed

	double tol[]       = { EPS, 2*EPS, EPS, 4*EPS, EPS, 2*EPS, EPS, };
	bool differences[] =
		{ diff_const_Plotting_Nodes(p_nodes_data_r->plt_line, p_nodes_data_c->plt_line, tol[0]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_tri,  p_nodes_data_c->plt_tri,  tol[1]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_quad, p_nodes_data_c->plt_quad, tol[2]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_tet,  p_nodes_data_c->plt_tet,  tol[3]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_hex,  p_nodes_data_c->plt_hex,  tol[4]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_wedge,p_nodes_data_c->plt_wedge,tol[5]),
		  diff_const_Plotting_Nodes(p_nodes_data_r->plt_pyr,  p_nodes_data_c->plt_pyr,  tol[6]),
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
		pass = false;

		if (differences[0])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_line, p_nodes_data_c->plt_line, tol[0]);
		if (differences[1])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_tri,  p_nodes_data_c->plt_tri,  tol[1]);
		if (differences[2])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_quad, p_nodes_data_c->plt_quad, tol[2]);
		if (differences[3])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_tet,  p_nodes_data_c->plt_tet,  tol[3]);
		if (differences[4])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_hex,  p_nodes_data_c->plt_hex,  tol[4]);
		if (differences[5])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_wedge,p_nodes_data_c->plt_wedge,tol[5]);
		if (differences[6])
			print_diff_const_Plotting_Nodes(p_nodes_data_r->plt_pyr,  p_nodes_data_c->plt_pyr,  tol[6]);
	}
	destructor_Plotting_Nodes_Data(p_nodes_data_r);
	destructor_Plotting_Nodes_Data(p_nodes_data_c);

	assert_condition(pass);
}

static void test_unit_nodes_face_correspondence (struct Test_Info*const test_info)
{
	UNUSED(test_info);
	char* test_name = "nodes_correspondence";
	bool pass = true;

	struct Nodes_FC_Data* nodes_fc_data_r = constructor_Nodes_FC_Data('r',test_name), // destructed
	                    * nodes_fc_data_c = constructor_Nodes_FC_Data('c',test_name); // destructed

	// Note: The TRI node ordering actually depends on the symmetries of the nodes (which is not being tested for
	//       here. Keep this in mind if additional triangular node sets are to be tested for face correspondence. This
	//       is not relevant for the tensor-product nodes.
	bool diffs[] =
		{ diff_const_Multiarray_Vector_i(nodes_fc_data_r->point_N,nodes_fc_data_c->point_N),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->line_3, nodes_fc_data_c->line_3),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->line_4, nodes_fc_data_c->line_4),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->quad_9, nodes_fc_data_c->quad_9),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->quad_16,nodes_fc_data_c->quad_16),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->tri_6,  nodes_fc_data_c->tri_6),
		  diff_const_Multiarray_Vector_i(nodes_fc_data_r->tri_10, nodes_fc_data_c->tri_10),
		};

	const int n_diff = sizeof(diffs)/sizeof(*diffs);
	if (check_diff(n_diff,diffs,&pass)) {
		int ind = 0;
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->point_N,nodes_fc_data_c->point_N);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->line_3, nodes_fc_data_c->line_3);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->line_4, nodes_fc_data_c->line_4);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->quad_9, nodes_fc_data_c->quad_9);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->quad_16,nodes_fc_data_c->quad_16);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->tri_6,  nodes_fc_data_c->tri_6);
		if (diffs[ind++]) print_diff_const_Multiarray_Vector_i(nodes_fc_data_r->tri_10, nodes_fc_data_c->tri_10);
	}
	destructor_Nodes_FC_Data(nodes_fc_data_r);
	destructor_Nodes_FC_Data(nodes_fc_data_c);

	assert_condition(pass);
}

// Level 1 ********************************************************************************************************** //

static struct Nodes_Data_TP* constructor_Nodes_Data_TP (const char eval_type, const char*const node_type)
{
	struct Nodes_Data_TP* nodes_data = calloc(1,sizeof *nodes_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"nodes/");
		strcat(file_name_part,node_type);
		const char*const file_name_full = set_data_file_name_unit(file_name_part);

		nodes_data->d1_p3_EQ  = constructor_file_name_const_Nodes("Nodes_d1_p3_EQ",file_name_full);   // keep
		nodes_data->d2_p4_EQ  = constructor_file_name_const_Nodes("Nodes_d2_p4_EQ",file_name_full);   // keep
		nodes_data->d3_p3_EQ  = constructor_file_name_const_Nodes("Nodes_d3_p3_EQ",file_name_full);   // keep
		nodes_data->d1_p3_GLL = constructor_file_name_const_Nodes("Nodes_d1_p3_GLoL",file_name_full); // keep
		nodes_data->d2_p4_GLL = constructor_file_name_const_Nodes("Nodes_d2_p4_GLoL",file_name_full); // keep
		nodes_data->d3_p3_GLL = constructor_file_name_const_Nodes("Nodes_d3_p3_GLoL",file_name_full); // keep
		nodes_data->d1_p3_GL  = constructor_file_name_const_Nodes("Nodes_d1_p3_GLe",file_name_full);  // keep
		nodes_data->d2_p4_GL  = constructor_file_name_const_Nodes("Nodes_d2_p4_GLe",file_name_full);  // keep
		nodes_data->d3_p3_GL  = constructor_file_name_const_Nodes("Nodes_d3_p3_GLe",file_name_full);  // keep
	} else if (eval_type == 'c') {
		nodes_data->d1_p3_EQ  = constructor_const_Nodes_tp(1,3,NODES_EQ);  // keep
		nodes_data->d2_p4_EQ  = constructor_const_Nodes_tp(2,4,NODES_EQ);  // keep
		nodes_data->d3_p3_EQ  = constructor_const_Nodes_tp(3,3,NODES_EQ);  // keep
		nodes_data->d1_p3_GLL = constructor_const_Nodes_tp(1,3,NODES_GLL); // keep
		nodes_data->d2_p4_GLL = constructor_const_Nodes_tp(2,4,NODES_GLL); // keep
		nodes_data->d3_p3_GLL = constructor_const_Nodes_tp(3,3,NODES_GLL); // keep
		nodes_data->d1_p3_GL  = constructor_const_Nodes_tp(1,3,NODES_GL);  // keep
		nodes_data->d2_p4_GL  = constructor_const_Nodes_tp(2,4,NODES_GL);  // keep
		nodes_data->d3_p3_GL  = constructor_const_Nodes_tp(3,3,NODES_GL);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return nodes_data;
}

static void destructor_Nodes_Data_TP (struct Nodes_Data_TP* nodes_data)
{
	destructor_const_Nodes(nodes_data->d1_p3_EQ);
	destructor_const_Nodes(nodes_data->d2_p4_EQ);
	destructor_const_Nodes(nodes_data->d3_p3_EQ);
	destructor_const_Nodes(nodes_data->d1_p3_GLL);
	destructor_const_Nodes(nodes_data->d2_p4_GLL);
	destructor_const_Nodes(nodes_data->d3_p3_GLL);
	destructor_const_Nodes(nodes_data->d1_p3_GL);
	destructor_const_Nodes(nodes_data->d2_p4_GL);
	destructor_const_Nodes(nodes_data->d3_p3_GL);

	free(nodes_data);
}

static struct Nodes_Data_SI* constructor_Nodes_Data_SI (const char eval_type, const char*const node_type)
{
	struct Nodes_Data_SI* nodes_data = calloc(1,sizeof *nodes_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"nodes/");
		strcat(file_name_part,node_type);
		const char*const file_name_full = set_data_file_name_unit(file_name_part);

		nodes_data->d2_p3_EQ  = constructor_file_name_const_Nodes("Nodes_d2_p3_EQ",file_name_full);  // keep
		nodes_data->d2_p3_AO  = constructor_file_name_const_Nodes("Nodes_d2_p3_AO",file_name_full);  // keep
		nodes_data->d2_p3_WSH = constructor_file_name_const_Nodes("Nodes_d2_p3_WSH",file_name_full); // keep
		nodes_data->d2_p6_WV  = constructor_file_name_const_Nodes("Nodes_d2_p6_WV",file_name_full);  // keep
		nodes_data->d3_p2_AO  = constructor_file_name_const_Nodes("Nodes_d3_p2_AO",file_name_full);  // keep
		nodes_data->d3_p2_WSH = constructor_file_name_const_Nodes("Nodes_d3_p2_WSH",file_name_full); // keep
		nodes_data->d3_p4_WV  = constructor_file_name_const_Nodes("Nodes_d3_p4_WV",file_name_full);  // keep
	} else if (eval_type == 'c') {
		nodes_data->d2_p3_EQ  = constructor_const_Nodes_si(2,3,NODES_EQ);  // keep
		nodes_data->d2_p3_AO  = constructor_const_Nodes_si(2,3,NODES_AO);  // keep
		nodes_data->d2_p3_WSH = constructor_const_Nodes_si(2,3,NODES_WSH); // keep
		nodes_data->d2_p6_WV  = constructor_const_Nodes_si(2,6,NODES_WV);  // keep
		nodes_data->d3_p2_AO  = constructor_const_Nodes_si(3,2,NODES_AO);  // keep
		nodes_data->d3_p2_WSH = constructor_const_Nodes_si(3,2,NODES_WSH); // keep
		nodes_data->d3_p4_WV  = constructor_const_Nodes_si(3,4,NODES_WV);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return nodes_data;
}

static void destructor_Nodes_Data_SI (struct Nodes_Data_SI* nodes_data)
{
	destructor_const_Nodes(nodes_data->d2_p3_EQ);
	destructor_const_Nodes(nodes_data->d2_p3_AO);
	destructor_const_Nodes(nodes_data->d2_p3_WSH);
	destructor_const_Nodes(nodes_data->d2_p6_WV);
	destructor_const_Nodes(nodes_data->d3_p2_AO);
	destructor_const_Nodes(nodes_data->d3_p2_WSH);
	destructor_const_Nodes(nodes_data->d3_p4_WV);

	free(nodes_data);
}

static struct Nodes_Data_PYR* constructor_Nodes_Data_PYR (const char eval_type, const char*const node_type)
{
	struct Nodes_Data_PYR* nodes_data = calloc(1,sizeof *nodes_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		strcpy(file_name_part,"nodes/");
		strcat(file_name_part,node_type);
		const char*const file_name_full = set_data_file_name_unit(file_name_part);

		nodes_data->d3_p2_GL  = constructor_file_name_const_Nodes("Nodes_d3_p2_GLe", file_name_full); // keep
		nodes_data->d3_p3_GLL = constructor_file_name_const_Nodes("Nodes_d3_p3_GLoL",file_name_full); // keep
		nodes_data->d3_p4_WV  = constructor_file_name_const_Nodes("Nodes_d3_p4_WV",  file_name_full); // keep
	} else if (eval_type == 'c') {
		nodes_data->d3_p2_GL  = constructor_const_Nodes_pyr(3,2,NODES_GL);  // keep
		nodes_data->d3_p3_GLL = constructor_const_Nodes_pyr(3,3,NODES_GLL); // keep
		nodes_data->d3_p4_WV  = constructor_const_Nodes_pyr(3,4,NODES_WV);  // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return nodes_data;
}

static void destructor_Nodes_Data_PYR (struct Nodes_Data_PYR* nodes_data)
{
	destructor_const_Nodes(nodes_data->d3_p2_GL);
	destructor_const_Nodes(nodes_data->d3_p3_GLL);
	destructor_const_Nodes(nodes_data->d3_p4_WV);

	free(nodes_data);
}

static struct Plotting_Nodes_Data* constructor_Plotting_Nodes_Data (const char eval_type, const char*const node_type)
{
	struct Plotting_Nodes_Data* p_nodes_data = calloc(1,sizeof *p_nodes_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		sprintf(file_name_part,"%s%s","nodes/",node_type);
		const char*const file_name_full = set_data_file_name_unit(file_name_part);

		p_nodes_data->plt_line  = constructor_file_name_const_Plotting_Nodes("Nodes_line", file_name_full); // keep
		p_nodes_data->plt_tri   = constructor_file_name_const_Plotting_Nodes("Nodes_tri",  file_name_full); // keep
		p_nodes_data->plt_quad  = constructor_file_name_const_Plotting_Nodes("Nodes_quad", file_name_full); // keep
		p_nodes_data->plt_tet   = constructor_file_name_const_Plotting_Nodes("Nodes_tet",  file_name_full); // keep
		p_nodes_data->plt_hex   = constructor_file_name_const_Plotting_Nodes("Nodes_hex",  file_name_full); // keep
		p_nodes_data->plt_wedge = constructor_file_name_const_Plotting_Nodes("Nodes_wedge",file_name_full); // keep
		p_nodes_data->plt_pyr   = constructor_file_name_const_Plotting_Nodes("Nodes_pyr",  file_name_full); // keep
	} else if (eval_type == 'c') {
		p_nodes_data->plt_line  = constructor_const_Plotting_Nodes(3,LINE);  // keep
		p_nodes_data->plt_tri   = constructor_const_Plotting_Nodes(3,TRI);   // keep
		p_nodes_data->plt_quad  = constructor_const_Plotting_Nodes(3,QUAD);  // keep
		p_nodes_data->plt_tet   = constructor_const_Plotting_Nodes(3,TET);   // keep
		p_nodes_data->plt_hex   = constructor_const_Plotting_Nodes(3,HEX);   // keep
		p_nodes_data->plt_wedge = constructor_const_Plotting_Nodes(3,WEDGE); // keep
		p_nodes_data->plt_pyr   = constructor_const_Plotting_Nodes(3,PYR);   // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return p_nodes_data;
}

static void destructor_Plotting_Nodes_Data (struct Plotting_Nodes_Data* p_nodes_data)
{
	destructor_const_Plotting_Nodes(p_nodes_data->plt_line);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_tri);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_quad);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_tet);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_hex);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_wedge);
	destructor_const_Plotting_Nodes(p_nodes_data->plt_pyr);

	free(p_nodes_data);
}

static struct Nodes_FC_Data* constructor_Nodes_FC_Data (const char eval_type, const char*const node_type)
{
	struct Nodes_FC_Data* nodes_fc_data = calloc(1,sizeof *nodes_fc_data); // returned
	if (eval_type == 'r') {
		char file_name_part[STRLEN_MAX];
		sprintf(file_name_part,"%s%s","nodes/",node_type);
		const char*const file_name_full = set_data_file_name_unit(file_name_part);

		nodes_fc_data->point_N = constructor_file_name_const_Multiarray_Vector_i("point_N",file_name_full); // keep
		nodes_fc_data->line_3  = constructor_file_name_const_Multiarray_Vector_i("line_3", file_name_full); // keep
		nodes_fc_data->line_4  = constructor_file_name_const_Multiarray_Vector_i("line_4", file_name_full); // keep
		nodes_fc_data->quad_9  = constructor_file_name_const_Multiarray_Vector_i("quad_9", file_name_full); // keep
		nodes_fc_data->quad_16 = constructor_file_name_const_Multiarray_Vector_i("quad_16",file_name_full); // keep
		nodes_fc_data->tri_6   = constructor_file_name_const_Multiarray_Vector_i("tri_6",  file_name_full); // keep
		nodes_fc_data->tri_10  = constructor_file_name_const_Multiarray_Vector_i("tri_10", file_name_full); // keep
	} else if (eval_type == 'c') {
		nodes_fc_data->point_N = constructor_nodes_face_corr(0,2,NODES_GL,ST_TP);  // keep
		nodes_fc_data->line_3  = constructor_nodes_face_corr(1,2,NODES_GL,ST_TP);  // keep
		nodes_fc_data->line_4  = constructor_nodes_face_corr(1,3,NODES_GL,ST_TP);  // keep
		nodes_fc_data->quad_9  = constructor_nodes_face_corr(2,2,NODES_GL,ST_TP);  // keep
		nodes_fc_data->quad_16 = constructor_nodes_face_corr(2,3,NODES_GL,ST_TP);  // keep
		nodes_fc_data->tri_6   = constructor_nodes_face_corr(2,2,NODES_WSH,ST_SI); // keep
		nodes_fc_data->tri_10  = constructor_nodes_face_corr(2,3,NODES_WSH,ST_SI); // keep
	} else {
		EXIT_UNSUPPORTED;
	}
	return nodes_fc_data;
}

static void destructor_Nodes_FC_Data (struct Nodes_FC_Data* nodes_fc_data)
{
	destructor_const_Multiarray_Vector_i(nodes_fc_data->point_N);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->line_3);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->line_4);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->quad_9);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->quad_16);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->tri_6);
	destructor_const_Multiarray_Vector_i(nodes_fc_data->tri_10);

	free(nodes_fc_data);
}

