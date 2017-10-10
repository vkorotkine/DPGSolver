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

#include "test_support_nodes.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support_matrix.h"
#include "test_support_multiarray.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_alloc.h"

#include "matrix.h"
#include "vector.h"

#include "file_processing.h"
#include "nodes.h"
#include "nodes_plotting.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Nodes\* from data in the input file of the given name.
 *  \return Standard. */
static struct Nodes* constructor_file_name_Nodes
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief Constructor for a \ref Plotting_Nodes container from data in the input file of the given name.
 *  \return Standard. */
struct Plotting_Nodes* constructor_file_name_Plotting_Nodes
	(const char*const var_name,      ///< Defined for \ref constructor_file_name_const_Plotting_Nodes.
	 const char*const file_name_full ///< Defined for \ref constructor_file_name_const_Plotting_Nodes.
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //

const struct const_Nodes* constructor_file_name_const_Nodes
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Nodes*) constructor_file_name_Nodes(var_name,file_name_full);
}

const struct const_Plotting_Nodes* constructor_file_name_const_Plotting_Nodes
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Plotting_Nodes*) constructor_file_name_Plotting_Nodes(var_name,file_name_full);
}

// Difference functions ********************************************************************************************* //

/** \brief Check the difference between members of the input \ref Nodes\*s.
 *  \return `true` if inputs differ; `false` otherwise. */
static bool diff_Nodes
	(const struct Nodes* a, ///< Input 0.
	 const struct Nodes* b, ///< Input 1.
	 const double tol       ///< The tolerance.
	);

/** \brief Check the difference between members of the input \ref Plotting_Nodes\*s.
 *  \return `true` if inputs differ; `false` otherwise. */
static bool diff_Plotting_Nodes
	(const struct Plotting_Nodes* a, ///< Input 0.
	 const struct Plotting_Nodes* b, ///< Input 1.
	 const double tol                ///< The tolerance.
	);

bool diff_const_Nodes (const struct const_Nodes*const a, const struct const_Nodes*const b, const double tol)
{
	return diff_Nodes((struct Nodes*)a,(struct Nodes*)b,tol);
}

bool diff_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const a, const struct const_Plotting_Nodes*const b, const double tol)
{
	return diff_Plotting_Nodes((const struct Plotting_Nodes*const)a,(const struct Plotting_Nodes*const)b,tol);
}

// Printing functions *********************************************************************************************** //

/// \brief Print the relative difference of the \ref Nodes members, outputting 0 if less than the tolerance.
static void print_diff_Nodes
	(const struct Nodes* a, ///< Input 0.
	 const struct Nodes* b, ///< Input 1.
	 const double tol       ///< The tolerance.
	);

/// \brief Print the relative difference of the \ref Plotting_Nodes members, outputting 0 if less than the tolerance.
static void print_diff_Plotting_Nodes
	(const struct Plotting_Nodes* a, ///< Input 0.
	 const struct Plotting_Nodes* b, ///< Input 1.
	 const double tol                ///< The tolerance.
	);

void print_diff_const_Nodes
	(const struct const_Nodes*const a, const struct const_Nodes*const b, const double tol)
{
	print_diff_Nodes((struct Nodes*)a,(struct Nodes*)b,tol);
}

void print_diff_const_Plotting_Nodes
	(const struct const_Plotting_Nodes*const a, const struct const_Plotting_Nodes*const b, const double tol)
{
	print_diff_Plotting_Nodes((struct Plotting_Nodes*)a,(struct Plotting_Nodes*)b,tol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for a \ref Nodes\* from the current line in the input file.
 *  \return Standard. */
static struct Nodes* constructor_file_Nodes
	(FILE* data_file ///< The pointer to the file from which to read the data.
	);

/** \brief Constructor for a \ref Plotting_Nodes container from the current line in the input file.
 *  \return Standard. */
static struct Plotting_Nodes* constructor_file_Plotting_Nodes
	(FILE* data_file ///< The pointer to the file from which to read the data.
	);

static struct Nodes* constructor_file_name_Nodes (const char*const var_name, const char*const file_name_full)
{
	struct Nodes* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Nodes(data_file);
		}
	}
	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

struct Plotting_Nodes* constructor_file_name_Plotting_Nodes (const char*const var_name, const char*const file_name_full)
{
	struct Plotting_Nodes* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Plotting_Nodes(data_file);
		}
	}
	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

static bool diff_Nodes (const struct Nodes* a, const struct Nodes* b, const double tol)
{
	if (diff_Matrix_d(a->rst,b->rst,tol))
		return true;

	if ((a->has_weights) && ((a->has_weights != b->has_weights) || diff_Vector_d(a->w,b->w,tol)))
		return true;

	return false;
}

static bool diff_Plotting_Nodes (const struct Plotting_Nodes* a, const struct Plotting_Nodes* b, const double tol)
{
	const struct Nodes* a_n = (const struct Nodes*)a,
	                  * b_n = (const struct Nodes*)b;

	if (diff_Nodes(a_n,b_n,tol))
		return true;

	if (diff_Multiarray_Vector_i(a->connect,b->connect))
		return true;

	if (diff_Multiarray_Vector_i(a->connect_e,b->connect_e))
		return true;

	if (diff_Vector_i(a->vtk_types,b->vtk_types))
		return true;

	return false;
}

static void print_diff_Nodes (const struct Nodes* a, const struct Nodes* b, const double tol)
{
	print_diff_Matrix_d(a->rst,b->rst,tol);

	printf("has weights: %d %d\n",a->has_weights,b->has_weights);
	if (a->has_weights && b->has_weights)
		print_diff_Vector_d(a->w,b->w,tol);
}

static void print_diff_Plotting_Nodes (const struct Plotting_Nodes* a, const struct Plotting_Nodes* b, const double tol)
{
	print_diff_Nodes((struct Nodes*)a,(struct Nodes*)b,tol);

	print_diff_Multiarray_Vector_i(a->connect,b->connect);
	print_diff_Multiarray_Vector_i(a->connect_e,b->connect_e);
	print_diff_Vector_i(a->vtk_types,b->vtk_types);
}

// Level 1 ********************************************************************************************************** //

static struct Nodes* constructor_file_Nodes (FILE* data_file)
{
	struct Nodes* nodes = calloc(1,sizeof *nodes); // returned

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	read_skip_file_const_b("has_weights",data_file,&nodes->has_weights);
	skip_lines(data_file,1);
	nodes->rst = constructor_file_Matrix_d(data_file,true); // keep

	if (nodes->has_weights) {
		skip_lines(data_file,1);
		struct Matrix_d* w_M = constructor_file_Matrix_d(data_file,true); // destructed
		nodes->w = constructor_move_Vector_d_Matrix_d(w_M); // keep
		destructor_Matrix_d(w_M);
	} else {
		nodes->w = NULL;
	}

	return nodes;
}

static struct Plotting_Nodes* constructor_file_Plotting_Nodes (FILE* data_file)
{
	struct Plotting_Nodes* p_nodes = calloc(1,sizeof *p_nodes); // returned
	struct Nodes* nodes = (struct Nodes*) p_nodes;

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	read_skip_file_const_b("has_weights",data_file,&nodes->has_weights);
	assert(!nodes->has_weights);
	skip_lines(data_file,1);
	nodes->rst = constructor_file_Matrix_d(data_file,true); // keep
	nodes->w   = NULL;

	skip_lines(data_file,1);
	p_nodes->connect = constructor_file_Multiarray_Vector_i(data_file); // keep

	skip_lines(data_file,1);
	p_nodes->connect_e = constructor_file_Multiarray_Vector_i(data_file); // keep

	skip_lines(data_file,1);
	p_nodes->vtk_types = constructor_file_Vector_i(data_file,true); // keep

	return p_nodes;
}
