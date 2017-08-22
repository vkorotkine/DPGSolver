// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/**	\file
 */

#include "file_processing.h"

#include <string.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>

#include "Macros.h"
#include "Parameters.h"

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

void skip_lines (FILE* file, char**const line, const int line_size, const int n_skip)
{
	for (int n = 0; n < n_skip; n++)
		fgets(*line,line_size,file);
}

void discard_line_values (char**const line, int n_discard)
{
	char* endptr = NULL;
	for (int n = 0; n < n_discard; n++) {
		strtod(*line,&endptr);
		*line = endptr;
	}
}

void read_line_values_i (char**const line, const ptrdiff_t n_val, int*const vals, const bool decrement)
{
	char* endptr = NULL;
	for (ptrdiff_t n = 0; n < n_val; n++) {
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
	for (ptrdiff_t n = 0; n < n_val; n++) {
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

void read_skip_i (const char*const line, int*const var)
{
	sscanf(line,"%*s %u",var);
}

void read_skip_const_c (const char*const line, const char*const var)
{
	sscanf(line,"%*s %s",(char*)var);
}

void read_skip_const_i (const char*const line, const int*const var)
{
	sscanf(line,"%*s %u",(int*)var);
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

void strcat_path_i (char* dest, const int src)
{
	char src_c[STRLEN_MIN] = {0};
	sprintf(src_c,"%u",src);
	strcat(dest,src_c);
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
