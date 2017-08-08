// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 *	\brief Brief file description if relevant.
 *
 *	General comments.
 */

#include "file_processing.h"

#include <string.h>

#include "Macros.h"
#include "Parameters.h"

FILE* fopen_checked (const char*const file_name_full)
{
	FILE* file = fopen(file_name_full,"r");
	if (file == NULL) {
		printf("File: '%s' is not present.\n",file_name_full);
		EXIT_UNSUPPORTED;
	}
	return file;
}

void discard_line_values (const char**const line, unsigned int n_discard)
{
	char* endptr = NULL;
	for (unsigned int n = 0; n < n_discard; n++) {
		strtod(*line,&endptr);
		*line = endptr;
	}
}

void read_line_values_ui (const char**const line, unsigned int n_val, unsigned int*const vals)
{
	char* endptr = NULL;
	for (unsigned int n = 0; n < n_val; n++) {
		vals[n] = strtol(*line,&endptr,10);
		*line = endptr;
	}
}

void read_skip_c (const char*const line, char*const var)
{
	sscanf(line,"%*s %s",var);
}

void read_skip_ui (const char*const line, unsigned int*const var)
{
	sscanf(line,"%*s %u",var);
}

void read_skip_const_c (const char*const line, const char*const var)
{
	sscanf(line,"%*s %s",(char*)var);
}

void read_skip_const_ui (const char*const line, const unsigned int*const var)
{
	sscanf(line,"%*s %u",(unsigned int*)var);
}

void read_skip_const_b (const char*const line, const bool*const var)
{
	sscanf(line,"%*s %d",(int*)var);
}

void strcat_path_c (char* dest, const char*const src, bool add_slash)
{
	if (!strstr(src,"NONE")) {
		strcat(dest,src);
		if (add_slash)
			strcat(dest,"/");
	}
}

void strcat_path_ui (char* dest, const unsigned int src)
{
	char src_c[STRLEN_MIN] = {0};
	sprintf(src_c,"%u",src);
	strcat(dest,src_c);
}
