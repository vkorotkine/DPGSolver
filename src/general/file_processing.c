// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "file_processing.h"

#include <string.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>

#include "macros.h"
#include "constants_alloc.h"

// Static function declarations ************************************************************************************* //

/** \brief Converts the string pointed to by the argument `str` to an `int`.
 *	\return See brief.
 *
 *	The implementation was taken from [this SO answer][SO_strtoi].
 *
 *	<!-- References: -->
 *	[SO_strtoi]: https://stackoverflow.com/a/29378380/5983549
 */
static int strtoi
	(const char* str, ///< The input string.
	 char** endptr,   ///< Will point to the next char following the `int` input in `str`.
	 int base         ///< The base for the number to be read (generally base 10).
	);

// Interface functions ********************************************************************************************** //

FILE* fopen_checked (const char*const file_name_full)
{
	FILE* file = fopen(file_name_full,"r");
	if (file == NULL) {
		printf("File: '%s' is not present.\n",file_name_full);
		EXIT_UNSUPPORTED;
	}
	return file;
}

FILE* fopen_input (const char*const input_path, const char*const input_spec)
{
	char input_name[STRLEN_MAX];

	strcpy(input_name,input_path);
	if (strstr(input_spec,"geometry"))
		strcat(input_name,"geometry_parameters.geo");
	else
		EXIT_UNSUPPORTED;

	return fopen_checked(input_name);
}

void skip_lines (FILE* file, const int n_skip)
{
	char line[STRLEN_MAX];
	for (int n = 0; n < n_skip; ++n)
		fgets(line,sizeof(line),file);
}

void skip_lines_ptr (FILE* file, char**const line, const int line_size, const int n_skip)
{
	for (int n = 0; n < n_skip; ++n)
		fgets(*line,line_size,file);
}

void discard_line_values (char**const line, int n_discard)
{
	char* endptr = NULL;
	for (int n = 0; n < n_discard; ++n) {
		strtod(*line,&endptr);
		*line = endptr;
	}
}

void read_line_values_i (char**const line, const ptrdiff_t n_val, int*const vals, const bool decrement)
{
	char* endptr = NULL;
	for (ptrdiff_t n = 0; n < n_val; ++n) {
		vals[n] = strtoi(*line,&endptr,10);
		*line = endptr;

		if (decrement) {
			if (!(vals[n] > 0))
				EXIT_UNSUPPORTED;
			--vals[n];
		}
	}
}

void read_line_values_l (char**const line, const ptrdiff_t n_val, long int*const vals, const bool decrement)
{
	char* endptr = NULL;
	for (ptrdiff_t n = 0; n < n_val; ++n) {
		vals[n] = strtol(*line,&endptr,10);
		*line = endptr;

		if (decrement) {
			if (!(vals[n] > 0))
				EXIT_UNSUPPORTED;
			--vals[n];
		}
	}
}

void read_line_values_d (char**const line, const ptrdiff_t n_val, double*const vals)
{
	char* endptr = NULL;
	for (ptrdiff_t n = 0; n < n_val; ++n) {
		vals[n] = strtod(*line,&endptr);
		*line = endptr;
	}
}

void read_skip_c (const char*const line, char*const var)
{
	sscanf(line,"%*s %s",var);
}

void read_skip_i (const char*const line, int*const var)
{
	sscanf(line,"%*s %u",var);
}

void read_skip_const_c (const char*const line, const char*const var)
{
	sscanf(line,"%*s %s",(char*)var);
}

void read_skip_const_b (const char*const line, const bool*const var)
{
	sscanf(line,"%*s %d",(int*)var);
}

void read_skip_const_d (char*const line, const double*const var, const int n_skip, const bool remove_semi)
{
	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	char* token_s_clean = ( remove_semi ? strtok(token_s,";") : token_s );

	sscanf(token_s_clean,"%lf",(double*)var);
}

void read_skip_const_i_1 (char*const line, const int n_skip, const int*const var, const int n_var)
{
	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip-1; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		token_s = strtok(NULL," ");
		sscanf(token_s,"%d",(int*)&var[i]);
	}
}

void read_skip_ptrdiff_1 (char*const line, const int n_skip, ptrdiff_t*const var, const int n_var)
{
	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip-1; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		token_s = strtok(NULL," ");
		sscanf(token_s,"%td",&var[i]);
	}
}

void read_skip_file_const_b (const char*const var_name, FILE* file, const bool*const var)
{
	char line[STRLEN_MAX];
	fgets(line,sizeof(line),file);

	if (!strstr(line,var_name))
		EXIT_ERROR("Did not find '%s' in the current line of the file.\n",var_name);

	read_skip_const_b(line,var);
}

void read_skip_file_const_i (const char*const var_name, FILE* file, const int*const var)
{
	char line[STRLEN_MAX];
	fgets(line,sizeof(line),file);

	if (!strstr(line,var_name))
		EXIT_ERROR("Did not find '%s' in the current line of the file.\n",var_name);

	read_skip_const_i_1(line,1,var,1);
}

void read_skip_file_i (const char*const var_name, FILE* file, int*const var)
{
	char line[STRLEN_MAX];
	fgets(line,sizeof(line),file);

	if (!strstr(line,var_name))
		EXIT_ERROR("Did not find '%s' in the current line (%s) of the file.\n",var_name,line);

	read_skip_i(line,var);
}

void strcat_path_c (char* dest, const char*const src, const char*const trail)
{
	if (!strstr(src,"NONE")) {
		strcat(dest,src);
		if (trail)
			strcat(dest,trail);
	}
}

void strcat_path_i (char* dest, const int src)
{
	char src_c[STRLEN_MIN] = {0};
	sprintf(src_c,"%u",src);
	strcat(dest,src_c);
}

char* extract_name (const char*const name_full, const bool extension_present)
{
	int len_name = 0;

	const char* beg = NULL, // Pointer to the first element of the name_full.
	          * end = NULL; // Pointer to one past the last element of the name_full.

	bool found_extension = false;
	int ind = 0;
	for (int i = strlen(name_full)-1; i >= 0; --i) {
		const char token = name_full[i];
		if (token == '/') {
			beg = &name_full[i+1];
			break;
		}

		if (!found_extension) {
			if (!extension_present) {
				end = &name_full[i+1];
				found_extension = true;
			} else if (token == '.') {
				end = &name_full[i];
				found_extension = true;
			}
		}

		if (!found_extension)
			continue;

		++len_name; // Note: One longer than the number of characters.
	}

	if (len_name == 0)
		EXIT_ERROR("Did not find the name in '%s'.\n",name_full);

	char*const name = calloc(len_name+1 , sizeof *name); // returned

	ind = 0;
	while (end != beg) {
		name[ind] = *beg++;
		++ind;
	}

	return name;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Converts the string pointed to by the argument `str` to a smaller type than `long int`.
 *	\return See brief.
 *
 *	The implementation was taken from [this SO answer][SO_strtoi].
 *
 *	<!-- References: -->
 *	[SO_strtoi]: https://stackoverflow.com/a/29378380/5983549
 */
static long strto_subrange
	(const char *str, ///< See strtoi.
	 char **endptr,   ///< See strtoi.
	 int base,        ///< See strtoi.
	 long min,        ///< The minimum supported value.
	 long max         ///< The maximum supported value.
	);

static int strtoi (const char* str, char** endptr, int base)
{
#if INT_MAX == LONG_MAX && INT_MIN == LONG_MIN
	return (int) strtol(str,endptr,base);
#else
	return (int) strto_subrange(str,endptr,base,INT_MIN,INT_MAX);
#endif
}

// Level 1 ********************************************************************************************************** //

static long strto_subrange (const char* str, char** endptr, int base, long min, long max)
{
	long y = strtol(str,endptr,base);
	if (y > max) {
		errno = ERANGE;
		return max;
	}
	if (y < min) {
		errno = ERANGE;
		return min;
	}
	return y;
}
