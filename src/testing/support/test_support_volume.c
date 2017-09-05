// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "test_support_volume.h"

#include "test_support_matrix.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "intrusive.h"
#include "volume.h"
#include "element.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Volume* constructor_Volume
	(FILE* file, const char*const line, const struct const_Intrusive_List*const elements)
{
	struct Volume* volume = malloc(sizeof *volume); // returned

	read_skip_const_i(line,&volume->index);
	read_skip_file_const_b("boundary",file,&volume->boundary);
	read_skip_file_const_b("curved",file,&volume->curved);
	skip_lines(file,1);
	struct Matrix_d* xyz_ve = constructor_file_Matrix_d(file,true); // keep
	const_constructor_move_Matrix_d(&volume->xyz_ve,xyz_ve);

	int elem_type;
	read_skip_file_i("elem_type",file,&elem_type);
	const_cast_const_Element(&volume->element,get_element_by_type(elements,elem_type));

	return volume;
}

struct Volume* get_volume_by_index (const struct Intrusive_List*const volumes, const int index)
{
	for (const struct Intrusive_Link* curr = volumes->first; curr; curr = curr->next) {
		struct Volume* volume = (struct Volume*) curr;
		if (volume->index == index)
			return volume;
	}
	EXIT_ERROR("Could not find the volume with index: %d.\n",index);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

