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

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_core.h"
#include "definitions_math.h"


#include "def_templates_geometry_parametric.h"

#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for the geometric data required for the supported parametric cases.
struct Geo_Data {
	Real x_scale,       ///< Multiplicative scaling for the x-coordinates.
	     r_i,           ///< Radius corresponding to the 'i'nner surface.
	     r_o,           ///< Radius corresponding to the 'o'uter surface.
	     total_radians, ///< The total number of radians to use when a cylindrical surface is present.
	     center[DMAX];  ///< The center of the geometry under consideration.
};

/** \brief Return the statically allocated \ref Geo_Data container.
 *  \return See brief. */
static struct Geo_Data get_geo_data
	(const char*const geo_name ///< The name corresponding to the type of geometry for which to obtain data.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_R* constructor_xyz_cylinder_parametric_T
	(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_R* xyz = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Real*const x_i = get_col_const_Multiarray_R(0,xyz_i),
	          *const y_i = get_col_const_Multiarray_R(1,xyz_i),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_R(2,xyz_i) : NULL );

	Real*const x = get_col_Multiarray_R(0,xyz),
	    *const y = get_col_Multiarray_R(1,xyz),
	    *const z = ( DIM > 2 ? get_col_Multiarray_R(2,xyz) : NULL );

	for (int n = 0; n < n_n; ++n) {
		const Real xy[2] = { x_i[n], y_i[n], };

		const Real r = GSL_MAX(fabs(xy[0]),fabs(xy[1]));

		Real t = atan2(xy[1],xy[0]);
		if (t >= -1.0*PI_OVER_4 && t < 1.0*PI_OVER_4)
			t = xy[1]/r*PI_OVER_4;
		else if ( t >=  1.0*PI_OVER_4 && t < 3.0*PI_OVER_4)
			t = 0.5*PI - xy[0]/r*PI_OVER_4;
		else if ((t >=  3.0*PI_OVER_4 && t <= 4.0*PI_OVER_4) || (t >= -4.0*PI_OVER_4 && t < -3.0*PI_OVER_4))
			t = PI - xy[1]/r*PI_OVER_4;
		else if (t >= -3.0*PI_OVER_4 && t < -1.0*PI_OVER_4)
			t = 1.5*PI + xy[0]/r*PI_OVER_4;
		else
			EXIT_ERROR("Unsupported: %f %f %f\n",xy[0],xy[1],t);

		x[n] = r*cos(t);
		y[n] = r*sin(t);
		if (DIM > 2)
			z[n] = z_i[n];
	}
	return (struct const_Multiarray_R*) xyz;
}

const struct const_Multiarray_R* constructor_xyz_trigonometric_cube_parametric_T
	(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_R* xyz = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Real*const x_i = get_col_const_Multiarray_R(0,xyz_i),
	          *const y_i = ( DIM > 1 ? get_col_const_Multiarray_R(1,xyz_i) : NULL ),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_R(2,xyz_i) : NULL );

	Real*const x = get_col_Multiarray_R(0,xyz),
	    *const y = ( DIM > 1 ? get_col_Multiarray_R(1,xyz) : NULL ),
	    *const z = ( DIM > 2 ? get_col_Multiarray_R(2,xyz) : NULL );

	const Real dxyz = 0.05;
	for (int n = 0; n < n_n; ++n) {
		if (DIM == 2) {
			x[n] = x_i[n] + dxyz*sin(PI*y_i[n]);
			y[n] = y_i[n];
		} else if (DIM == 3) {
			EXIT_ADD_SUPPORT; UNUSED(z_i); UNUSED(z);
		} else {
			EXIT_UNSUPPORTED;
		}
	}
	return (struct const_Multiarray_R*) xyz;
}

const struct const_Multiarray_R* constructor_xyz_joukowski_parametric_T
	(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(n_type);
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const ptrdiff_t n_n = xyz_i->extents[0];

	struct Multiarray_R* xyz = constructor_empty_Multiarray_R('C',2,(ptrdiff_t[]){n_n,DIM}); // returned

	const Real*const x_i = get_col_const_Multiarray_R(0,xyz_i),
	          *const y_i = get_col_const_Multiarray_R(1,xyz_i),
	          *const z_i = ( DIM > 2 ? get_col_const_Multiarray_R(2,xyz_i) : NULL );

	Real*const x = get_col_Multiarray_R(0,xyz),
	    *const y = get_col_Multiarray_R(1,xyz),
	    *const z = ( DIM > 2 ? get_col_Multiarray_R(2,xyz) : NULL );

	struct Geo_Data geo_data = get_geo_data("joukowski");

	const Real x_scale       = geo_data.x_scale,
	           r_i           = geo_data.r_i,
	           r_o           = geo_data.r_o,
	           total_radians = geo_data.total_radians,
		    *const center  = geo_data.center;

	for (int n = 0; n < n_n; ++n) {
		const Real t = 0.5*total_radians*(1.0-x_i[n]);

		const Complex zeta = center[0]+r_i*cos(t) + I*(center[1]+r_i*sin(t));
		const Complex z_j = zeta + 1.0/zeta;

		const Real x_j = creal(z_j)*x_scale,
		           y_j = cimag(z_j);

		const Real x_c = r_o*cos(t) - x_scale,
		           y_c = r_o*sin(t);

		x[n] = x_j*(1-y_i[n]) + x_c*(y_i[n]);
		y[n] = y_j*(1-y_i[n]) + y_c*(y_i[n]);
		if (DIM > 2)
			z[n] = z_i[n];
	}
	return (struct const_Multiarray_R*) xyz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Read the required geometry data for the Joukowski parametric domain into the \ref Geo_Data.
static void read_data_joukowski
	(struct Geo_Data*const geo_data ///< \ref Geo_Data.
	);

static struct Geo_Data get_geo_data (const char*const geo_name)
{
	static bool need_input = true;
	static struct Geo_Data geo_data;
	if (need_input) {
		need_input = false;
		if (strcmp(geo_name,"joukowski") == 0)
			read_data_joukowski(&geo_data);
		else
			EXIT_ERROR("Unsupported: %s.\n",geo_name);
	}

	return geo_data;
}

// Level 1 ********************************************************************************************************** //

static void read_data_joukowski (struct Geo_Data*const geo_data)
{
	const int count_to_find = 5;

	FILE* input_file = fopen_input('g',NULL,NULL); // closed

	Real total_degrees = 0.0;

	int count_found = 0;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_d("r_i", &count_found,line,&geo_data->r_i);
		read_skip_string_count_d("r_o", &count_found,line,&geo_data->r_o);
		read_skip_string_count_d("x_scale",&count_found,line,&geo_data->x_scale);
		read_skip_string_count_d("total_degrees",&count_found,line,&total_degrees);
		if (strstr(line,"center_cyl")) {
			read_skip_d_1(line,1,geo_data->center,DMAX);
			++count_found;
		}
	}
	fclose(input_file);

	geo_data->total_radians = total_degrees*PI/180.0;
	assert(count_found == count_to_find);
}
