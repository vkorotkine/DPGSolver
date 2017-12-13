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

#include "test_support.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"

// Static function declarations ************************************************************************************* //



// Interface functions ********************************************************************************************** //

void update_pass (bool*const pass, const bool pass_new)
{
	*pass = ( *pass == false ? *pass : pass_new );
}

bool compare_i (const int a, const int b, const bool print_enabled, const char*const var_name)
{
	bool pass = true;

	if (a != b) {
		pass = false;
		if (print_enabled)
			printf("%s: %d %d\n",var_name,a,b);
	}
	return pass;
}

bool compare_b (const bool a, const bool b, const bool print_enabled, const char*const var_name)
{
	bool pass = true;

	if (a != b) {
		pass = false;
		if (print_enabled)
			printf("%s: %d %d\n",var_name,a,b);
	}
	return pass;
}

char* set_print_name_container_member (const char*const name_container, int ind_container, const char*const name_member)
{
	static char print_name[STRLEN_MAX];

	strcpy(print_name,name_container);
	strcat(print_name,"[");

	char ind_str[STRLEN_MIN];
	sprintf(ind_str,"%d",ind_container);
	strcat(print_name,ind_str);

	strcat(print_name,"]::");
	strcat(print_name,name_member);

	return print_name;
}

void check_container_type (FILE* data_file, const char*const container_type)
{
	char line[STRLEN_MAX];
	if (fgets(line,sizeof(line),data_file) != NULL) {};

	char expected_line[STRLEN_MAX];
	strcpy(expected_line,"container ");
	strcat(expected_line,container_type);

	const bool found = ( strstr(line,expected_line) ? true : false );
	if (!found)
		EXIT_ERROR("Reading incorrect container type: %s. (expected: %s)",line,expected_line);
}

bool check_diff (const int n_entries, const bool*const differences, bool* pass)
{
	*pass = true;
	for (int i = 0; i < n_entries; ++i) {
		if (differences[i]) {
			*pass = false;
			printf("%d\n",i);
		}
	}
	return !(*pass);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

