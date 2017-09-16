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

#include "test_support_cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_alloc.h"

#include "matrix.h"
#include "vector.h"

#include "cubature.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Cubature\* from the current line in the input file.
 *  \return Standard. */
static struct Cubature* constructor_file_Cubature
	(FILE* data_file ///< The pointer to the file from which to read the data.
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //

struct Cubature* constructor_file_name_Cubature (const char*const var_name, const char*const file_name_full)
{
	struct Cubature* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Cubature(data_file);
		}
	}
	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

const struct const_Cubature* constructor_file_name_const_Cubature
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Cubature*) constructor_file_name_Cubature(var_name,file_name_full);
}

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between members of the input \ref Cubature\*s.
 *  \return The `true` if inputs differ; `false` otherwise. */
static bool diff_Cubature
	(const struct Cubature* a, ///< Input 0.
	 const struct Cubature* b, ///< Input 1.
	 const double tol          ///< The tolerance.
	);

bool diff_const_Cubature (const struct const_Cubature*const a, const struct const_Cubature*const b, const double tol)
{
	return diff_Cubature((struct Cubature*)a,(struct Cubature*)b,tol);
}

// Printing functions *********************************************************************************************** //

/// \brief Print the relative difference of the \ref Cubature members, outputting 0 if less than the tolerance.
static void print_diff_Cubature
	(const struct Cubature* a, ///< Input 0.
	 const struct Cubature* b, ///< Input 1.
	 const double tol          ///< The tolerance.
	);

void print_diff_const_Cubature
	(const struct const_Cubature*const a, const struct const_Cubature*const b, const double tol)
{
	print_diff_Cubature((struct Cubature*)a,(struct Cubature*)b,tol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Cubature* constructor_file_Cubature (FILE* data_file)
{
	struct Cubature* cubature = calloc(1,sizeof *cubature); // returned

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	read_skip_file_const_b("has_weights",data_file,&cubature->has_weights);
	skip_lines(data_file,1);
	cubature->rst = constructor_file_Matrix_d(data_file,true); // keep

	if (cubature->has_weights) {
		skip_lines(data_file,1);
		struct Matrix_d* w_M = constructor_file_Matrix_d(data_file,true); // destructed
		cubature->w = constructor_move_Vector_d_Matrix_d(w_M); // keep
		destructor_Matrix_d(w_M);
	} else {
		cubature->w = NULL;
	}

	return cubature;
}

static bool diff_Cubature (const struct Cubature* a, const struct Cubature* b, const double tol)
{
	if (diff_Matrix_d(a->rst,b->rst,tol))
		return true;

	if ((a->has_weights) && ((a->has_weights != b->has_weights) || diff_Vector_d(a->w,b->w,tol)))
		return true;

	return false;
}

static void print_diff_Cubature (const struct Cubature* a, const struct Cubature* b, const double tol)
{
	print_diff_Matrix_d(a->rst,b->rst,tol);

	printf("has weights: %d %d\n",a->has_weights,b->has_weights);
	if (a->has_weights && b->has_weights)
		print_diff_Vector_d(a->w,b->w,tol);
}
