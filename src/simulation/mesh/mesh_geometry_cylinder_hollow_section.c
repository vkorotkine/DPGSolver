// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "mesh_geometry_cylinder_hollow_section.h"
#include "vector.h"
#include "matrix.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Macros.h"
#include "constants_mesh.h"
#include "constants_alloc.h"

#include "math_functions.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for geometry data relating to 'c'ylindrical 'h'ollow 's'ections.
struct Geom_Data__chs {
	const double r_i, ///< The inner radius.
	             r_o; ///< The outer radius.
};

static struct Geom_Data__chs read_data_cylinder__hollow_section
	(const char*const input_path ///< Defined in \ref fopen_input.
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
		geom_data  = read_data_cylinder__hollow_section(input_path);
	}

	// Snap vertices to the boundary
	const int dims_to_check = 2;

	const ptrdiff_t n_n = nodes->extents[0];
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

static struct Geom_Data__chs read_data_cylinder__hollow_section (const char*const input_path)
{
	struct Geom_Data__chs geom_data;

	int       count_found   = 0;
	const int count_to_find = 2;

	FILE* input_file = fopen_input(input_path,"geometry"); // closed

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		if (strstr(line,"r_i")) {
			++count_found;
			read_skip_const_d(line,&geom_data.r_i,2,true);
		} else if (strstr(line,"r_o")) {
			++count_found;
			read_skip_const_d(line,&geom_data.r_o,2,true);
		}
	}
	fclose(input_file);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");

	return geom_data;
}
