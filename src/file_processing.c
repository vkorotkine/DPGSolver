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

void skip_lines (FILE* file, char**const line, const size_t line_size, const unsigned int n_skip)
{
	for (unsigned int n = 0; n < n_skip; n++)
		fgets(*line,line_size,file);
}

void discard_line_values (char**const line, unsigned int n_discard)
{
	char* endptr = NULL;
	for (unsigned int n = 0; n < n_discard; n++) {
		strtod(*line,&endptr);
		*line = endptr;
	}
}

void read_line_values_ui (char**const line, const unsigned int n_val, unsigned int*const vals, const bool decrement)
{
	char* endptr = NULL;
	for (unsigned int n = 0; n < n_val; n++) {
		vals[n] = strtol(*line,&endptr,10);
		*line = endptr;

		if (decrement) {
			if (!(vals[n] > 0))
				EXIT_UNSUPPORTED;
			--vals[n];
		}
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
