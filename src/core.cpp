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

#include "core.h"

#include <assert.h>

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_core.h"

#include "file_processing.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

const char* set_petsc_options_name (const char* petsc_options_name)
{
	static char full_name[STRLEN_MAX] = { 0, };
	sprintf(full_name,"%s%s%s%s",PROJECT_SOURCE_DIR,"external/petsc/options_files/",petsc_options_name,".txt");

	FILE* file = fopen_checked(full_name);
	fclose(file);

	return full_name;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
