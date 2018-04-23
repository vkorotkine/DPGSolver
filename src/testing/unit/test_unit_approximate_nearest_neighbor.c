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

#include "test_base.h"
#include "test_support.h"
#include "test_support_matrix.h"

#include "macros.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for 'A'pproximate 'N'earest 'N'eighbor info.
struct ANN_Info {
	// Read from file.
	const struct const_Matrix_d* nodes_b;  ///< Background nodes.
	const struct const_Matrix_d* nodes_s;  ///< Search nodes.
	/** Indices of approximate nearest neighbors in \ref ANN_Info::nodes_b for each of the nodes in
	 *  \ref ANN_Info::nodes_s. */
	const struct const_Vector_i* inds_ann;
};

/** \brief Constructor for a \ref ANN_Info container.
 *  \return See brief. */
static const struct ANN_Info* constructor_ANN_Info
	(const char*const data_name ///< String containing the name of the data file (without path).
	);

/// \brief Destructor for a \ref ANN_Info container.
static void destructor_ANN_Info
	(struct ANN_Info* ann_info ///< Standard.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs unit testing for the approximate nearest neighbor functionality
 *        (\ref test_unit_approximate_nearest_neighbor.c).
 *  \return 0 on success.
 *
 *  The test checks that the nearest neighbor algorithm finds the expected nearest neighbors (whose indices are
 *  specified in the input file).
 *
 *  \todo Optionally remove the expected indices and use a naive O(n^2) algorithm to find the nearest points.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	assert_condition_message(argc == 2,"Invalid number of input arguments");

	const struct ANN_Info*const ann_info = constructor_ANN_Info(argv[1]); // destructed

print_const_Matrix_d(ann_info->nodes_b);
EXIT_ADD_SUPPORT;
	destructor_ANN_Info((struct ANN_Info*)ann_info);

	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct ANN_Info* constructor_ANN_Info (const char*const data_name)
{
	char d_name[STRLEN_MAX];
	sprintf(d_name,"%s%s","ann/",data_name);
	const char*const f_name = set_data_file_name_unit(d_name);

	struct ANN_Info*const ann_info = calloc(1,sizeof(*ann_info)); // free

	ann_info->nodes_b  = constructor_file_name_const_Matrix_d("nodes_background",f_name); // destructed
	ann_info->nodes_s  = constructor_file_name_const_Matrix_d("nodes_search",    f_name); // destructed
	ann_info->inds_ann = constructor_file_name_const_Vector_i("indices_ann",     f_name); // destructed

	return ann_info;
}

static void destructor_ANN_Info (struct ANN_Info* ann_info)
{
	destructor_const_Matrix_d(ann_info->nodes_b);
	destructor_const_Matrix_d(ann_info->nodes_s);
	destructor_const_Vector_i(ann_info->inds_ann);
	free(ann_info);
}
