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

#include "test_support_volume.h"

#include "test_support_multiarray.h"

#include <stdlib.h>
#include <stdio.h>

#include "macros.h"

#include "multiarray.h"
#include "matrix.h"

#include "intrusive.h"
#include "volume.h"
#include "element.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

struct Volume* constructor_Volume
	(FILE* file, char* line, const struct const_Intrusive_List*const elements)
{
	struct Volume* volume = calloc(1,sizeof *volume); // returned

	read_skip_const_i_1(line,1,&volume->index,1);
	read_skip_file_const_b("boundary",file,&volume->boundary);
	read_skip_file_const_b("curved",file,&volume->curved);
	skip_lines(file,1);
	struct Multiarray_d* xyz_ve = constructor_file_Multiarray_d(file,true); // keep
	const_constructor_move_Multiarray_d(&volume->xyz_ve,xyz_ve);

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

