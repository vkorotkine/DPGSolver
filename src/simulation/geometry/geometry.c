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

#include "geometry.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "computational_elements.h"
#include "volume_solver.h"
#include "face_solver.h"

#include "const_cast.h"
#include "element_solver.h"
#include "file_processing.h"
#include "geometry_normals.h"
#include "intrusive.h"
#include "operator.h"
#include "multiarray_operator.h"
#include "simulation.h"
#include "test_case.h"
#include "visualization.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "geometry_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "geometry_T.c"
#include "undef_templates_type.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

bool is_internal_geom_straight ( )
{
	static bool igs = false;
	static bool need_input = true;
	if (need_input) {
		char line[STRLEN_MAX];
		FILE* input_file = fopen_input('t',NULL,NULL); // closed
		while (fgets(line,sizeof(line),input_file)) {
			if (strstr(line,"use_straight_internal_geometry"))
				read_skip_const_b(line,&igs);
		}
		fclose(input_file);
	}
	return igs;
}

bool geometry_depends_on_face_pointers()
{
	if (get_set_domain_type(NULL) && is_internal_geom_straight())
		return true;
	return false;
}
