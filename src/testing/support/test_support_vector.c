// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support.h"
#include "test_support_multiarray.h"

#include "macros.h"
#include "definitions_alloc.h"

#include "file_processing.h"
#include "vector.h"
#include "math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //
// Constructor functions ******************************************************************************************** //

struct Vector_d* constructor_file_name_Vector_d (const char*const var_name, const char*const file_name_full)
{
	struct Vector_d* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file) != NULL) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Vector_d(data_file,true);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

struct Vector_i* constructor_file_name_Vector_i (const char*const var_name, const char*const file_name_full)
{
	struct Vector_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file) != NULL) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_file_Vector_i(data_file,true);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
}

const struct const_Vector_d* constructor_file_name_const_Vector_d
	(const char*const var_name, const char*const file_name_full)
{
	return (const struct const_Vector_d*) constructor_file_name_Vector_d(var_name,file_name_full);
}


struct Vector_d* constructor_file_Vector_d (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Vector_d");

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	ptrdiff_t ext_0 = 0;
	read_line_values_l(line_ptr,1,&ext_0,false);

	double data[ext_0];
	read_line_values_d(line_ptr,ext_0,data);

	return constructor_copy_Vector_d_d(ext_0,data);
}

struct Vector_i* constructor_file_Vector_i (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Vector_i");

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	ptrdiff_t ext_0 = 0;
	read_line_values_l(line_ptr,1,&ext_0,false);

	int data[ext_0];
	read_line_values_i(line_ptr,ext_0,data,false);

	return constructor_copy_Vector_i_i(ext_0,data);
}

// Difference functions ********************************************************************************************* //

bool diff_Vector_d (const struct Vector_d*const a, const struct Vector_d*const b, const double tol)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0)
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (!equal_d(a->data[i],b->data[i],tol))
			return true;
	}
	return false;
}

bool diff_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0)
		return true;

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (a->data[i] != b->data[i])
			return true;
	}
	return false;
}

// Printing functions *********************************************************************************************** //

void print_diff_Vector_d (const struct Vector_d*const a, const struct Vector_d*const b, const double tol)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0) {
		printf("Note: Attempting to compare Vectors of different size:\n");
		printf("a (%td):\n",a->ext_0);
		print_Vector_d(a,tol);
		printf("b (%td):\n",b->ext_0);
		print_Vector_d(b,tol);
		return;
	}

	const double*const data_a = a->data,
	            *const data_b = b->data;

	double data[size];
	for (ptrdiff_t i = 0; i < size; ++i)
		data[i] = norm_diff_d(1,&data_a[i],&data_b[i],"Inf");

	// Temporarily modify the data pointer to print the vector of differences without constructing a new object.
	double* data_ptr = a->data;

	struct Vector_d* a_tmp = (struct Vector_d*) a;
	a_tmp->data = data;
	print_Vector_d(a_tmp,tol);
	a_tmp->data = data_ptr;
}

void print_diff_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0) {
		printf("Note: Attempting to compare Vectors of different size:\n");
		print_Vector_i(a);
		print_Vector_i(b);
		return;
	}

	const int*const data_a = a->data,
	         *const data_b = b->data;

	for (ptrdiff_t i = 0; i < size; i++) {
		printf("% 12d ",data_a[i]-data_b[i]);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

