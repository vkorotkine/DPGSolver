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

#include "geometry_normals.h"

#include "multiarray.h"

#include "face_solver.h"

#include "file_processing.h"
#include "geometry_parametric.h"
#include "math_functions.h"
#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "geometry_normals_T.c"

bool using_exact_normals ( )
{
	static bool need_input = true;
	static bool flag       = false;
	if (need_input) {
		need_input = false;
		char line[STRLEN_MAX];
		FILE* input_file = input_file = fopen_input('t',NULL,NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"use_exact_normals_for_all")) read_skip_const_b(line,&flag);
		}
		fclose(input_file);
	}
	return flag;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
