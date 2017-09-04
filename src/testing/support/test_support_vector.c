// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support_multiarray.h"

#include "Macros.h"
#include "file_processing.h"
#include "vector.h"

#include "constants_alloc.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Vector_i* constructor_file_name_Vector_i (const char*const var_name, const char*const file_name_full)
{
	struct Vector_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
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

struct Vector_i* constructor_file_Vector_i (FILE* data_file, const bool check_container)
{
	if (check_container)
		check_container_type(data_file,"Vector_i");

	char line[STRLEN_MAX];
	char* line_ptr[1] = {line};
	fgets(line,sizeof(line),data_file);

	ptrdiff_t ext_0 = 0;
	read_line_values_l(line_ptr,1,&ext_0,false);

	int data[ext_0];
	read_line_values_i(line_ptr,ext_0,data,false);

	return constructor_copy_Vector_i_i(ext_0,data);
}

bool diff_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;

	if (size != b->ext_0)
		EXIT_ERROR("Comparing Vectors of different size");

	for (ptrdiff_t i = 0; i < size; ++i) {
		if (a->data[i] != b->data[i])
			return true;
	}

	return false;
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

