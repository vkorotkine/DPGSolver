// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_multiarray.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Macros.h"
#include "file_processing.h"

#include "constants_alloc.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref Multiarray_Vector_i as read from a file.
 *	\return Standard. */
static struct Multiarray_Vector_i* constructor_read_Multiarray_Vector_i
	(FILE* data_file ///< The file containing the data.
	);

// Interface functions ********************************************************************************************** //

struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i
	(const char*const var_name, const char*const file_name_full)
{
	struct Multiarray_Vector_i* dest = NULL;

	FILE* data_file = fopen_checked(file_name_full); // closed

	bool found_var = false;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),data_file)) {
		if (strstr(line,var_name)) {
			found_var = true;
			dest = constructor_read_Multiarray_Vector_i(data_file);
		}
	}

	fclose(data_file);

	if (!found_var)
		EXIT_ERROR("Did not find the '%s' variable in the file: %s",var_name,file_name_full);

	return dest;
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

/// \brief Check that the container type is that which is expected.
static void check_container_type
	(FILE* data_file,                ///< The file containing the data.
	 const char*const container_type ///< The container type.
	);

/** \brief Obtain the order and extents of the Multiarray.
 *	\return See brief. */
struct Multiarray_Partial read_order_extents
	(FILE* data_file ///< The file containing the data.
	);

static struct Multiarray_Vector_i* constructor_read_Multiarray_Vector_i (FILE* data_file)
{

	check_container_type(data_file,"Multiarray_Vector_i");

	struct Multiarray_Partial ma_p = read_order_extents(data_file);

printf("%d %td\n",ma_p.order,ma_p.extents[0]);
EXIT_UNSUPPORTED;

/// \todo Continue here.
	struct Multiarray_Vector_i* dest = NULL;

	return dest;
}

// Level 1 ********************************************************************************************************** //

static void check_container_type (FILE* data_file, const char*const container_type)
{
	char line[STRLEN_MAX];
	fgets(line,sizeof(line),data_file);

	if (strstr(container_type,"Multiarray_Vector_i")) {
		printf("%s\n",line);
		if (!strstr(line,"container Multiarray_Vector_i"))
			EXIT_ERROR("Reading incorrect container type: %s",line);
	} else {
		EXIT_UNSUPPORTED;
	}
}

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
