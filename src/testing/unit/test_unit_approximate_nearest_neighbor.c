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

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the approximate nearest neighbor functionality
 *        (\ref test_unit_approximate_nearest_neighbor.c).
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

	const struct Input_ANN*const ann_i    = constructor_file_Input_ANN(argv[1]);  // destructed
	const struct const_Vector_i* nn_naive = constructor_ann_indices_naive(ann_i); // destructed
	const struct const_Vector_i* nn_nlogn = constructor_ann_indices(ann_i);       // destructed

	bool pass = true;
	if (diff_const_Vector_i(nn_naive,nn_nlogn)) {
		pass = false;

		print_const_Matrix_d(ann_i->nodes_b);
		print_const_Matrix_d(ann_i->nodes_s);
		print_diff_const_Vector_i(nn_naive,nn_nlogn);
		print_const_Vector_i(nn_naive);
		print_const_Vector_i(nn_nlogn);
	}
	assert_condition(pass);

	destructor_Input_ANN((struct Input_ANN*)ann_i);
	destructor_const_Vector_i(nn_naive);
	destructor_const_Vector_i(nn_nlogn);

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
