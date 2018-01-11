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

#include "file_processing_conversions.h"

#include <assert.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_geometry.h"

#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Get the value of the macro defined constant associated with the string variable input of the definition type.
 *  \return See brief. */
static int get_define
	(const char*const def_str, ///< String associated with macro variable definition.
	 const char*const def_type ///< Defined for \ref read_skip_convert_const_i.
	);

// Interface functions ********************************************************************************************** //

void read_skip_convert_const_i
	(const char*const line, const char*const def_type, const int*const var, int*const count_found)
{
	if (!strstr(line,def_type))
		return;

	if (count_found)
		++*count_found;

	char def_str[STRLEN_MAX];
	sscanf(line,"%*s %s",def_str);

	const int var_def = get_define(def_str,def_type);
	const_cast_i(var,var_def);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static int get_define (const char*const def_str, const char*const def_type)
{
	int def_i = -1;
	if (strcmp(def_type,"geom_parametrization") == 0) {
		if      (strcmp(def_str,"radial_proj") == 0) def_i = GEOM_PRM_RADIAL_PROJ;
		else if (strcmp(def_str,"arc_length")  == 0) def_i = GEOM_PRM_ARC_LENGTH;
		else if (strcmp(def_str,"normal_proj") == 0) def_i = GEOM_PRM_NORMAL_PROJ;
		else                                         EXIT_ERROR("Unsupported: %s\n",def_str);
	} else {
		EXIT_ERROR("Unsupported: %s\n",def_type);
	}

	return def_i;
}
