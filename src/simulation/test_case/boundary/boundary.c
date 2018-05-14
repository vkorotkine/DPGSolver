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

#include "boundary.h"

#include <string.h>
#include "macros.h"

#include "definitions_bc.h"

#include "multiarray.h"
#include "vector.h"

#include "element_solver.h"
#include "face.h"
#include "face_solver.h"
#include "volume.h"
#include "volume_solver.h"

#include "compute_face_rlhs.h"
#include "file_processing.h"
#include "geometry_normals.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solve_dg.h"

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "boundary_T.c"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

bool using_exact_normals_for_boundary ( )
{
	static bool need_input = true;
	static bool flag       = false;
	if (need_input) {
		need_input = false;
		char line[STRLEN_MAX];
		FILE* input_file = input_file = fopen_input('t',NULL,NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"use_exact_normals_for_boundary")) read_skip_const_b(line,&flag);
		}
		fclose(input_file);

		if (flag && using_exact_normals()) {
			printf("*** Warning: Enabling exact boundary normals is redundant if exact normals are enabled. *** \n");
			PRINT_FILELINE;
		}
	}
	return flag;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
