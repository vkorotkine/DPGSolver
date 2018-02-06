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
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_core.h"
#include "definitions_math.h"


#include "def_templates_geometry_parametric.h"

#include "def_templates_volume_solver.h"

#include "def_templates_multiarray.h"

// Static function declarations ************************************************************************************* //

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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
