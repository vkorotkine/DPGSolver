// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
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
	fgets(line,sizeof(line),data_file);

	char expected_line[STRLEN_MAX];
	strcpy(expected_line,"container ");
	strcat(expected_line,container_type);

	const bool found = ( strstr(line,expected_line) ? true : false );
	if (!found)
		EXIT_ERROR("Reading incorrect container type: %s",line);
}

bool check_diff (const int n_entries, const bool*const differences)
{
	for (int i = 0; i < n_entries; ++i) {
		if (differences[i])
			return true;
	}
	return false;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

