// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_multiarray.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "test_support_vector.h"

#include "Macros.h"
#include "file_processing.h"
#include "multiarray.h"

#include "constants_alloc.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Multiarray_Vector_i\* as read from a file.
 *	\return Standard. */
static struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i
	(FILE* data_file ///< The file containing the data.
	);

// Interface functions ********************************************************************************************** //

struct Multiarray_Vector_i* constructor_file_name_Multiarray_Vector_i
	(const char*const var_name, const char*const file_name_full)
{
	struct Multiarray_Vector_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
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

void check_container_type (FILE* data_file, const char*const container_type)
{
	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);

	char expected_line[STRLEN_MAX];
	strcpy(expected_line,"container ");
	strcat(expected_line,container_type);

	const bool found = ( strstr(line,expected_line) ? true : false );
	if (!found)
		EXIT_ERROR("Reading incorrect container type: %s",line);
}

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

	switch (order) {
	case 1:
		for (ptrdiff_t i = 0; i < size; i++)
			print_diff_Vector_i(a->data[i],b->data[i]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
	printf("\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

static struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i (FILE* data_file)
{
	check_container_type(data_file,"Multiarray_Vector_i");

	struct Multiarray_Partial ma_p = read_order_extents(data_file);

	struct Multiarray_Vector_i* dest = NULL;
	switch (ma_p.order) {
	case 1:
		dest = constructor_empty_Multiarray_Vector_i(false,ma_p.order,ma_p.extents[0]);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);
	if (!strstr(line,"ext_0/data"))
		EXIT_ERROR("Did not find expected data description.");

	ptrdiff_t size = compute_size(dest->order,dest->extents);
	for (ptrdiff_t i = 0; i < size; ++i)
		dest->data[i] = constructor_file_Vector_i(data_file,false);

	return dest;
}

// Level 1 ********************************************************************************************************** //

struct Multiarray_Partial read_order_extents (FILE* data_file)
{
	struct Multiarray_Partial ma_p;

	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);
	read_skip_i(line,&ma_p.order);

	if (ma_p.order > EXTENTS_MAX)
		EXIT_ERROR("Increase the size of EXTENTS_MAX or use dynamically allocated extents.");

	fgets(line,sizeof(line),data_file);
	read_skip_ptrdiff_1(line,1,ma_p.extents,ma_p.order);

	return ma_p;
}
