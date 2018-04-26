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
/// \file

#include "restart.h"

#include <string.h>

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "volume_solver.h"

#include "approximate_nearest_neighbor.h"
#include "computational_elements.h"
#include "const_cast.h"
#include "element.h"
#include "file_processing.h"
#include "intrusive.h"
#include "simulation.h"
#include "solution.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "restart_T.c"

const char* get_restart_name ( )
{
	static char restart_path[STRLEN_MAX];
	static bool needs_input = true;

	if (needs_input) {
		needs_input = false;
		FILE* input_file = fopen_input('c',NULL,NULL); // closed
		char line[STRLEN_MAX];
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"restart_path"))
				read_skip_c(line,restart_path);
		}
		fclose(input_file);
	}

	return restart_path;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
