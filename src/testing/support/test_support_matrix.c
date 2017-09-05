// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_matrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support_multiarray.h"

#include "Macros.h"
#include "file_processing.h"
#include "math_functions.h"
#include "matrix.h"

#include "constants_alloc.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Matrix_i* constructor_file_name_Matrix_i (const char*const var_name, const char*const file_name_full)
{
	struct Matrix_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Matrix_i(data_file,true);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s.",var_name,file_name_full);

	return dest;
}
struct Matrix_d* constructor_file_name_Matrix_d (const char*const var_name, const char*const file_name_full)
{
	struct Matrix_d* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Matrix_d(data_file,true);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s.",var_name,file_name_full);

	return dest;
}

struct Matrix_i* constructor_file_Matrix_i (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Matrix_i");

	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);

	char layout = 0;
	ptrdiff_t ext_0 = 0,
	          ext_1 = 0;

	sscanf(line,"%c %td %td",&layout,&ext_0,&ext_1);

	const ptrdiff_t size = ext_0*ext_1;

	// Read data by row and transpose if necessary.
	int data[size];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		fgets(line,sizeof(line),data_file);

		char* line_ptr[1] = {line};
		read_line_values_i(line_ptr,ext_1,&data[i*ext_1],false);
	}

	struct Matrix_i* dest = constructor_copy_Matrix_i_i('R',ext_0,ext_1,data);

	if (layout == 'C')
		EXIT_ADD_SUPPORT; // Implement matrix transpose (with fixed extents).

	return dest;
}

struct Matrix_d* constructor_file_Matrix_d (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Matrix_d");

	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);

	char layout = 0;
	ptrdiff_t ext_0 = 0,
	          ext_1 = 0;

	sscanf(line,"%c %td %td",&layout,&ext_0,&ext_1);

	const ptrdiff_t size = ext_0*ext_1;

	// Read data by row and transpose if necessary.
	double data[size];
	for (ptrdiff_t i = 0; i < ext_0; ++i) {
		fgets(line,sizeof(line),data_file);

		char* line_ptr[1] = {line};
		read_line_values_d(line_ptr,ext_1,&data[i*ext_1]);
	}

	struct Matrix_d* dest = constructor_copy_Matrix_d_d('R',ext_0,ext_1,data);

	if (layout == 'C')
		EXIT_ADD_SUPPORT; // Implement matrix transpose (with fixed extents).

	return dest;
}

bool diff_Matrix_i (const struct Matrix_i*const a, const struct Matrix_i*const b)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);

	if ((size != (b->ext_0)*(b->ext_1)) || (a->layout != b->layout))
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (a->data[i] != b->data[i])
			return true;
	}
	return false;
}

bool diff_Matrix_d (const struct Matrix_d*const a, const struct Matrix_d*const b, const double tol)
{
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);

	if ((size != (b->ext_0)*(b->ext_1)) || (a->layout != b->layout))
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (!equal_d(a->data[i],b->data[i],tol))
			return true;
	}
	return false;
}

bool diff_const_Matrix_d (const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, const double tol)
{
	return diff_Matrix_d((const struct Matrix_d*const)a,(const struct Matrix_d*const)b,tol);
}

void print_diff_Matrix_i (const struct Matrix_i*const a, const struct Matrix_i*const b)
{
	const char layout    = a->layout;
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);

	if ((size != (b->ext_0)*(b->ext_1)) || (layout != b->layout)) {
		printf("Attempting to compare Matrices of different size:\n");
		print_Matrix_i(a);
		print_Matrix_i(b);
		return;
	}

	const int*const data_a = a->data,
	         *const data_b = b->data;

	int data[size];
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] =  data_a[i]-data_b[i];

	// Temporarily modify the data pointer to print the matrix of differences without constructing a new object.
	int* data_ptr = a->data;

	struct Matrix_i* a_tmp = (struct Matrix_i*) a;
	a_tmp->data = data;
	print_Matrix_i(a_tmp);
	a_tmp->data = data_ptr;
}

void print_diff_Matrix_d (const struct Matrix_d*const a, const struct Matrix_d*const b, const double tol)
{
	const char layout    = a->layout;
	const ptrdiff_t size = (a->ext_0)*(a->ext_1);

//	if (1||(size != (b->ext_0)*(b->ext_1)) || (layout != b->layout)) {
	if ((size != (b->ext_0)*(b->ext_1)) || (layout != b->layout)) {
		printf("Attempting to compare Matrices of different size:\n");
		print_Matrix_d(a,tol);
		print_Matrix_d(b,tol);
		return;
	}

	const double*const data_a = a->data,
	            *const data_b = b->data;

	double data[size];
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = norm_diff_d(1,&data_a[i],&data_b[i],"Inf");

	// Temporarily modify the data pointer to print the matrix of differences without constructing a new object.
	double* data_ptr = a->data;

	struct Matrix_d* a_tmp = (struct Matrix_d*) a;
	a_tmp->data = data;
	print_Matrix_d(a_tmp,tol);
	a_tmp->data = data_ptr;
}

void print_diff_const_Matrix_d
	(const struct const_Matrix_d*const a, const struct const_Matrix_d*const b, const double tol)
{
	print_diff_Matrix_d((const struct Matrix_d*const)a,(const struct Matrix_d*const)b,tol);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

