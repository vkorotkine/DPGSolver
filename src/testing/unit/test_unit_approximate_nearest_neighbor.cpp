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

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "vector.h"

#include "approximate_nearest_neighbor.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for 'A'pproximate 'N'earest 'N'eighbor info.
struct ANN_Info {
	// Read from file.
	const struct const_Matrix_d* nodes_b;  ///< Background nodes.
	const struct const_Matrix_d* nodes_s;  ///< Search nodes.
};

/** \brief Constructor for a \ref Input_ANN container.
 *  \return See brief. */
static const struct Input_ANN* constructor_file_Input_ANN
	(const char*const data_name ///< String containing the name of the data file (without path).
	);

/** \brief Constructor for the approximate nearest neighbor indices using the naive O(n^2) algorithm.
 *  \return See brief. */
static const struct const_Vector_i* constructor_ann_indices_naive
	(const struct Input_ANN*const ann_i ///< Standard.
	);

/** \brief Constructor for a \ref Nodes_Sorted_ANN container using the naive O(n^2) nearest neighbor algorithm.
 *  \return See brief. */
static const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_naive
	(const struct const_Matrix_d*const nodes_i ///< Input unsorted nodes.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the approximate nearest neighbor functionality
 *        (\ref test_unit_approximate_nearest_neighbor.cpp).
 *  \return 0 on success.
 *
 *  The test checks that the nearest neighbor algorithm finds the expected nearest neighbors as determined by the naive
 *  O(n^2) algorithm.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const struct Input_ANN*const ann_i = constructor_file_Input_ANN(argv[1]); // destructed
	const struct const_Vector_i*const nn_naive = constructor_ann_indices_naive(ann_i), // destructed
	                           *const nn_ann   = constructor_ann_indices(ann_i);       // destructed
	const struct Nodes_Sorted_ANN*const nodes_sorted_naive = constructor_Nodes_Sorted_naive(ann_i->nodes_b), // dest.
	                             *const nodes_sorted_ann   = constructor_Nodes_Sorted_ANN(ann_i->nodes_b);   // dest.

	bool pass = true;
	const double tol = EPS;
	const bool* differences = (bool[])
		{ diff_const_Vector_i(nn_naive,nn_ann),
		  diff_const_Vector_i(nodes_sorted_naive->indices,nodes_sorted_ann->indices),
		  diff_const_Matrix_d(nodes_sorted_naive->nodes,nodes_sorted_ann->nodes,tol),
		};
	if (check_diff(3,differences,&pass)) {
		print_const_Matrix_d(ann_i->nodes_b);
		print_const_Matrix_d(ann_i->nodes_s);
		if (differences[0]) {
			print_diff_const_Vector_i(nn_naive,nn_ann);
			print_const_Vector_i(nn_naive);
			print_const_Vector_i(nn_ann);
		}
		if (differences[1]) {
			print_diff_const_Vector_i(nodes_sorted_naive->indices,nodes_sorted_ann->indices);
			print_const_Vector_i(nodes_sorted_naive->indices);
			print_const_Vector_i(nodes_sorted_ann->indices);
		}
		if (differences[2]) {
			print_diff_const_Matrix_d(nodes_sorted_naive->nodes,nodes_sorted_ann->nodes,tol);
			print_const_Matrix_d(nodes_sorted_naive->nodes);
			print_const_Matrix_d(nodes_sorted_ann->nodes);
		}
	}
	assert_condition(pass);

	destructor_Input_ANN((struct Input_ANN*)ann_i);
	destructor_const_Vector_i(nn_naive);
	destructor_const_Vector_i(nn_ann);
	destructor_Nodes_Sorted_ANN(nodes_sorted_naive);
	destructor_Nodes_Sorted_ANN(nodes_sorted_ann);

	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct Input_ANN* constructor_file_Input_ANN (const char*const data_name)
{
	char d_name[STRLEN_MAX];
	sprintf(d_name,"%s%s","ann/",data_name);
	const char*const f_name = set_data_file_name_unit(d_name);

	struct Input_ANN*const ann_i = calloc(1,sizeof(*ann_i)); // free

	ann_i->nodes_b = constructor_file_name_const_Matrix_d("nodes_background",f_name); // destructed
	ann_i->nodes_s = constructor_file_name_const_Matrix_d("nodes_search",    f_name); // destructed

	return ann_i;
}

static const struct const_Vector_i* constructor_ann_indices_naive (const struct Input_ANN*const ann_i)
{
	const struct const_Matrix_d*const nodes_b = ann_i->nodes_b;
	const struct const_Matrix_d*const nodes_s = ann_i->nodes_s;
	assert(nodes_b->ext_1 == DIM);
	assert(nodes_s->ext_1 == DIM);

	const ptrdiff_t n_b = nodes_b->ext_0,
	                n_s = nodes_s->ext_0;

	struct Matrix_d*const d2 = constructor_zero_Matrix_d('R',n_s,n_b); // destructed
	struct Vector_i*const nn = constructor_zero_Vector_i(n_s);         // returned
	for (int i = 0; i < n_s; ++i) {
		const double*const data_s = get_row_const_Matrix_d(i,nodes_s);

		double*const data_d2 = get_row_Matrix_d(i,d2);
		for (int j = 0; j < n_b; ++j) {
			const double*const data_b = get_row_const_Matrix_d(j,nodes_b);
			for (int d = 0; d < DIM; ++d)
				data_d2[j] += POW2(data_s[d]-data_b[d]);

			if (data_d2[j] < data_d2[nn->data[i]])
				nn->data[i] = j;
		}
	}
	destructor_Matrix_d(d2);

	return (struct const_Vector_i*) nn;
}

static const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_naive (const struct const_Matrix_d*const nodes_i)
{
	const ptrdiff_t n_n = nodes_i->ext_0;
	struct Vector_i*const inds_sorted  = constructor_empty_Vector_i(n_n); // moved
	int*const data_is = inds_sorted->data;
	for (int i = 0; i < n_n; ++i)
		inds_sorted->data[i] = i;

	struct Matrix_d nodes_s = { .layout = 'R', .ext_0 = 1, .ext_1 = DIM, .owns_data = false, .data = NULL, };
	struct Matrix_d*const nodes_b = constructor_copy_Matrix_d((struct Matrix_d*)nodes_i); // moved
	double*const data_b0 = nodes_b->data;

	assert(nodes_b->ext_1 == DIM);
	for (int ind_s = 1; ind_s < n_n; ++ind_s) {
		nodes_s.data   = get_row_Matrix_d(0,nodes_b);
		nodes_b->data  = get_row_Matrix_d(1,nodes_b);
		nodes_b->ext_0 -= 1;
		inds_sorted->data  += 1;
		inds_sorted->ext_0 -= 1;

		struct Input_ANN ann_i = { .nodes_b = (struct const_Matrix_d*) nodes_b,
		                           .nodes_s = (struct const_Matrix_d*) &nodes_s, };
		const struct const_Vector_i* ind_n = constructor_ann_indices_naive(&ann_i); // destructed

		swap_rows_Matrix_d(nodes_b,0,ind_n->data[0]);
		swap_vals_Vector_i(inds_sorted,0,ind_n->data[0]);
		destructor_const_Vector_i(ind_n);
	}

	nodes_b->ext_0 = n_n;
	nodes_b->data  = data_b0;

	inds_sorted->ext_0 = n_n;
	inds_sorted->data  = data_is;

	struct Nodes_Sorted_ANN*const nsa = calloc(1,sizeof *nsa); // returned
	nsa->nodes   = (struct const_Matrix_d*) nodes_b;     // keep
	nsa->indices = (struct const_Vector_i*) inds_sorted; // keep

	return nsa;
}
