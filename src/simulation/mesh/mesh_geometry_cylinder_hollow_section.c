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

#include "mesh_geometry_cylinder_hollow_section.h"
#include "vector.h"
#include "matrix.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_mesh.h"
#include "definitions_alloc.h"

#include "math_functions.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for geometry data relating to 'c'ylindrical 'h'ollow 's'ections.
struct Geom_Data__chs {
	const double r_i, ///< The inner radius.
	             r_o; ///< The outer radius.
};

/// \brief Read the required geometry data into \ref Geom_Data__chs.
static void read_data_cylinder__hollow_section
	(const char*const input_path,          ///< Defined in \ref fopen_input.
	 struct Geom_Data__chs*const geom_data ///< \ref Geom_Data__chs.
	);

// Interface functions ********************************************************************************************** //

void mesh_snap_to_cylinder__hollow_section
	(const char*const input_path, const struct const_Vector_i*const ve_curved, const struct Matrix_d*const nodes)
{
	// Set geometry data
	static bool need_input = true;

	static struct Geom_Data__chs geom_data;
	if (need_input) {
		need_input = false;
		read_data_cylinder__hollow_section(input_path,&geom_data);
	}

	// Snap vertices to the boundary
	const int dims_to_check = 2;

	const ptrdiff_t n_n = nodes->ext_0;
	for (ptrdiff_t n = 0; n < n_n; ++n) {
		if (!ve_curved->data[n])
			continue;

		const double r = norm_d(dims_to_check,get_row_Matrix_d(n,nodes),"L2");

		double r_ex = 0.0;
		if (equal_d(r,geom_data.r_i,NODETOL_MESH))
			r_ex = geom_data.r_i;
		else if (equal_d(r,geom_data.r_o,NODETOL_MESH))
			r_ex = geom_data.r_o;

		if (r_ex == 0.0)
			EXIT_ERROR("Did not find a matching curved surface");

		double*const ve_xyz = get_row_Matrix_d(n,nodes);
		const double t = atan2(ve_xyz[1],ve_xyz[0]);

		ve_xyz[0] = r_ex*cos(t);
		ve_xyz[1] = r_ex*sin(t);
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void read_data_cylinder__hollow_section (const char*const input_path, struct Geom_Data__chs*const geom_data)
{
	int       count_found   = 0;
	const int count_to_find = 2;

	FILE* input_file = fopen_input(input_path,'g'); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"r_i")) {
			++count_found;
			read_skip_const_d(line,&geom_data->r_i,2,true);
		} else if (strstr(line,"r_o")) {
			++count_found;
			read_skip_const_d(line,&geom_data->r_o,2,true);
		}
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}
