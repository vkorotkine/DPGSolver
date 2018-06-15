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

#include "file_processing.h"
#include "multiarray.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "geometry_parametric_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "geometry_parametric_T.c"
#include "undef_templates_type.h"

double f_gaussian_bump (const double x, const int diff_degree, const struct Function_Data_GP*const f_data)
{
	const struct Geo_Data geo_data = get_geo_data("gaussian_bump");

	const double scale = f_data->scale;

	const double a = geo_data.a,
	             b = geo_data.b,
	             c = geo_data.c,
	             d = geo_data.d;

	switch (diff_degree) {
	case 0:
		return scale*(1+d*(x-b))*a*exp(-1.0*pow(x-b,2.0)/(2.0*c*c));
		break;
	case 1: {
		const double df0 = (1+d*(x-b)),
		             df1 = d,
		             dg0 = a*exp(-1.0*pow(x-b,2.0)/(2.0*c*c)),
		             dg1 = a*(b-x)/(c*c)*exp(-0.5*pow((x-b)/c,2.0));
		return scale*(df1*dg0+df0*dg1);
		break;
	} default:
		EXIT_ERROR("Add support.\n");
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
