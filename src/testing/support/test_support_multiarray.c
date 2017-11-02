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

#include "test_support_multiarray.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support.h"
#include "test_support_matrix.h"
#include "test_support_vector.h"

#include "macros.h"
#include "definitions_alloc.h"

#include "multiarray.h"

#include "file_processing.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

///\{ \name The maximum number of extents for the read Multiarray.
#define EXTENTS_MAX 10
///\}

/// Container for partial members of the Multiarray.
struct Multiarray_Partial {
	int       order;                ///< Defined in \ref Multiarray_d.
	ptrdiff_t extents[EXTENTS_MAX]; ///< Defined in \ref Multiarray_d.
};

/** \brief Obtain the order and extents of the Multiarray.
 *	\return See brief. */
struct Multiarray_Partial read_order_extents
	(FILE* data_file ///< The file containing the data.
	);

/** \brief Constructor for a \ref Multiarray_Matrix_d\* as read from a file.
 *  \return Standard. */
static struct Multiarray_Matrix_d* constructor_file_Multiarray_Matrix_d
	(FILE* data_file ///< The file containing the data.
	);

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //

struct Multiarray_d* constructor_file_name_Multiarray_d
	(const char*const var_name, const char*const file_name_full)
{
	struct Multiarray_d* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file) != NULL) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Multiarray_d(data_file,true);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

const struct const_Multiarray_d* constructor_file_name_const_Multiarray_d
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Multiarray_d*) constructor_file_name_Multiarray_d(var_name,file_name_full);
}

struct Multiarray_Vector_i* constructor_file_name_Multiarray_Vector_i
	(const char*const var_name, const char*const file_name_full)
{
	struct Multiarray_Vector_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file) != NULL) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Multiarray_Vector_i(data_file);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

const struct const_Multiarray_Vector_i* constructor_file_name_const_Multiarray_Vector_i
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Multiarray_Vector_i*)
		constructor_file_name_Multiarray_Vector_i(var_name,file_name_full);
}

const struct const_Multiarray_Matrix_d* constructor_file_name_const_Multiarray_Matrix_d
	(const char*const var_name, const char*const file_name_full)
{
	struct Multiarray_Matrix_d* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file) != NULL) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Multiarray_Matrix_d(data_file);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return (const struct const_Multiarray_Matrix_d*) dest;
}

struct Multiarray_d* constructor_file_Multiarray_d (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Multiarray_d");

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	char layout = 0;
	int  order  = 0;
	sscanf(line,"%c %d",&layout,&order);

	if (order == 0)
		EXIT_UNSUPPORTED;

	ptrdiff_t extents[order];
	read_skip_ptrdiff_1 (line,2,extents,order);

	const ptrdiff_t size = compute_size(order,extents);

	// Read data by row and transpose if necessary.
	if (order > 2)
		EXIT_ADD_SUPPORT;

	ptrdiff_t ext_0 = extents[0],
	          ext_1 = extents[1];

	double* data = malloc(size * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		if (fgets(line,sizeof(line),data_file) != NULL) {};

		char* line_ptr[1] = {line};
		read_line_values_d(line_ptr,ext_1,&data[i*ext_1]);
	}

	struct Multiarray_d* dest = constructor_move_Multiarray_d_d('R',order,extents,true,data); // returned

	if (layout == 'C')
		transpose_Multiarray_d(dest,true);

	return dest;
}

struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i (FILE* data_file)
{
	check_container_type(data_file,"Multiarray_Vector_i");

	struct Multiarray_Partial ma_p = read_order_extents(data_file);

	struct Multiarray_Vector_i* dest = NULL;
	switch (ma_p.order) {
	case 1:
		dest = constructor_empty_Multiarray_Vector_i(false,ma_p.order,ma_p.extents);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};
	if (!strstr(line,"ext_0/data"))
		EXIT_ERROR("Did not find expected data description.");

	ptrdiff_t size = compute_size(dest->order,dest->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		dest->data[i] = constructor_file_Vector_i(data_file,false);

	return dest;
}

// Difference functions ********************************************************************************************* //

bool diff_Multiarray_Vector_i (const struct Multiarray_Vector_i*const a, const struct Multiarray_Vector_i*const b)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if (size != compute_size(b->order,b->extents))
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (diff_Vector_i(a->data[i],b->data[i]))
			return true;
	}

	return false;
}

bool diff_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a, const struct const_Multiarray_Vector_i*const b)
{
	return diff_Multiarray_Vector_i((struct Multiarray_Vector_i*)a,(struct Multiarray_Vector_i*)b);
}

bool diff_Multiarray_d (const struct Multiarray_d*const a, const struct Multiarray_d*const b, const double tol)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if ((size != compute_size(b->order,b->extents)) || (a->layout != b->layout))
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (!equal_d(a->data[i],b->data[i],tol))
			return true;
	}
	return false;
}

bool diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, const struct const_Multiarray_d*const b, const double tol)
{
	return diff_Multiarray_d((struct Multiarray_d*)a,(struct Multiarray_d*)b,tol);
}

bool diff_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, const struct Multiarray_Matrix_d*const b, const double tol)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if (size != compute_size(b->order,b->extents))
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (diff_Matrix_d(a->data[i],b->data[i],tol))
			return true;
	}

	return false;
}

bool diff_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, const struct const_Multiarray_Matrix_d*const b,
	 const double tol)
{
	return diff_Multiarray_Matrix_d((struct Multiarray_Matrix_d*)a,(struct Multiarray_Matrix_d*)b,tol);
}

// Printing functions *********************************************************************************************** //

void print_diff_Multiarray_Vector_i (const struct Multiarray_Vector_i*const a, const struct Multiarray_Vector_i*const b)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if (size != compute_size(b->order,b->extents)) {
		printf("Note: Attempting to compare Multiarrays of different size:\n");
		print_Multiarray_Vector_i(a);
		print_Multiarray_Vector_i(b);
		return;
	}

	const int order                = a->order;
	const ptrdiff_t *const extents = a->extents;

	printf("(diff) Multi-array extents: {");
	for (ptrdiff_t i = 0; i < order; i++)
		printf(" %zu,",extents[i]);
	printf(" }\n\n");

	for (ptrdiff_t i = 0; i < size; i++)
		print_diff_Vector_i(a->data[i],b->data[i]);
	printf("\n");
}

void print_diff_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a, const struct const_Multiarray_Vector_i*const b)
{
	print_diff_Multiarray_Vector_i((struct Multiarray_Vector_i*)a,(struct Multiarray_Vector_i*)b);
}

void print_diff_Multiarray_d (const struct Multiarray_d*const a, const struct Multiarray_d*const b, const double tol)
{
	const char layout    = a->layout;
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if ((size != compute_size(b->order,b->extents)) || (layout != b->layout)) {
		printf("Attempting to compare Multiarrays of different size:\n");
		print_Multiarray_d_tol(a,tol);
		print_Multiarray_d_tol(b,tol);
		return;
	}

	const double*const data_a = a->data,
	            *const data_b = b->data;

	double data[size];
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = norm_diff_d(1,&data_a[i],&data_b[i],"Inf");

	// Temporarily modify the data pointer to print the matrix of differences without constructing a new object.
	double* data_ptr = a->data;

	struct Multiarray_d* a_tmp = (struct Multiarray_d*) a;
	a_tmp->data = data;
	print_Multiarray_d_tol(a_tmp,tol);
	a_tmp->data = data_ptr;
}

void print_diff_const_Multiarray_d
	(const struct const_Multiarray_d*const a, const struct const_Multiarray_d*const b, const double tol)
{
	print_diff_Multiarray_d((const struct Multiarray_d*const)a,(const struct Multiarray_d*const)b,tol);
}

void print_diff_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, const struct Multiarray_Matrix_d*const b, const double tol)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);

	if (size != compute_size(b->order,b->extents)) {
		printf("Note: Attempting to compare Multiarrays of different size:\n");
		print_Multiarray_Matrix_d_tol(a,tol);
		print_Multiarray_Matrix_d_tol(b,tol);
		return;
	}

	const int order                = a->order;
	const ptrdiff_t *const extents = a->extents;

	printf("(diff) Multi-array extents: {");
	for (ptrdiff_t i = 0; i < order; i++)
		printf(" %zu,",extents[i]);
	printf(" }\n\n");

	for (ptrdiff_t i = 0; i < size; i++)
		print_diff_Matrix_d(a->data[i],b->data[i],tol);
	printf("\n");
}

void print_diff_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, const struct const_Multiarray_Matrix_d*const b,
	 const double tol)
{
	print_diff_Multiarray_Matrix_d((struct Multiarray_Matrix_d*)a,(struct Multiarray_Matrix_d*)b,tol);
}

// Math functions *************************************************************************************************** //

void perturb_Multiarray_d (struct Multiarray_d* a, const double da)
{
	const ptrdiff_t size = compute_size(a->order,a->extents);
	for (int i = 0; i < size; ++i)
		a->data[i] += da*(((double) rand())/((double) RAND_MAX+1));
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Multiarray_Matrix_d* constructor_file_Multiarray_Matrix_d (FILE* data_file)
{
	check_container_type(data_file,"Multiarray_Matrix_d");

	struct Multiarray_Partial ma_p = read_order_extents(data_file);

	struct Multiarray_Matrix_d* dest = NULL;
	switch (ma_p.order) {
	case 1:
		dest = constructor_empty_Multiarray_Matrix_d(false,ma_p.order,ma_p.extents);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	ptrdiff_t size = compute_size(dest->order,dest->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		dest->data[i] = constructor_file_Matrix_d(data_file,true);

	return dest;
}

struct Multiarray_Partial read_order_extents (FILE* data_file)
{
	struct Multiarray_Partial ma_p;

	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};
	read_skip_i(line,&ma_p.order);

	if (ma_p.order > EXTENTS_MAX)
		EXIT_ERROR("Increase the size of EXTENTS_MAX or use dynamically allocated extents.");

	if (fgets(line,sizeof(line),data_file) != NULL) {};
	read_skip_ptrdiff_1(line,1,ma_p.extents,ma_p.order);

	return ma_p;
}
