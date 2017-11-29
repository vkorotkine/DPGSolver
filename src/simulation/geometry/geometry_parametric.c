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

#include "geometry_parametric.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_math.h"

#include "multiarray.h"

// Static function declarations ************************************************************************************* //

#define PI_OVER_4 PI/4.0 ///< Pi divided by 4.

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_d* constructor_xyz_cylinder_parametric
	(const struct const_Multiarray_d* xyz_i, const struct Solver_Volume* s_vol, const struct Simulation* sim)
{
	UNUSED(s_vol);
	UNUSED(sim);
	const ptrdiff_t n_n = xyz_i->extents[0],
	                d   = xyz_i->extents[1];
	assert(d >= 2);

	struct Multiarray_d* xyz = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_n,d}); // returned

	const double*const x_i = get_col_const_Multiarray_d(0,xyz_i),
	            *const y_i = get_col_const_Multiarray_d(1,xyz_i),
	            *const z_i = ( d > 2 ? get_col_const_Multiarray_d(2,xyz_i) : NULL );

	double*const x = get_col_Multiarray_d(0,xyz),
	      *const y = get_col_Multiarray_d(1,xyz),
	      *const z = ( d > 2 ? get_col_Multiarray_d(2,xyz) : NULL );

	for (int n = 0; n < n_n; ++n) {
		const double xy[2] = { x_i[n], y_i[n], };

		const double r = GSL_MAX(xy[0],xy[1]);

		double t = atan2(xy[1],xy[0]);
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
		if (z)
			z[n] = z_i[n];
	}
	return (struct const_Multiarray_d*) xyz;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
