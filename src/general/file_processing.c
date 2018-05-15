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

#include "file_processing.h"

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

#include "macros.h"
#include "definitions_alloc.h"

#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/** \brief Converts the string pointed to by the argument `str` to an `int`.
 *  \return See brief.
 *
 *  The implementation was taken from [this SO answer][SO_strtoi].
 *
 *  <!-- References: -->
 *  [SO_strtoi]: https://stackoverflow.com/a/29378380/5983549
 */
static int strtoi
	(const char* str, ///< The input string.
	 char** endptr,   ///< Will point to the next char following the `int` input in `str`.
	 int base         ///< The base for the number to be read (generally base 10).
	);

/** \brief Analogue of `mkdir -p` but called from within the code.
 *  \return 0 if successful.
 *
 *  This function was taken from [Jonathon Reinhart's improved version] of [this SO answer][SO_mkdir_p].
 *
 *  <!-- References: -->
 *  [SO_mkdir_p]: https://stackoverflow.com/a/2336245/5983549
 *  [Jonathon Reinhart's improved version]: https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950
 */
static int mkdir_p
	(const char *path ///< The full path to the directory to be created.
	);

/** \brief Return a statically allocated `char*` holding the full name of the 's'erial/'p'arallel file.
 *  \return See brief. */
static const char* get_file_name_sp
	(const char sp_type,              ///< Defined for \ref fopen_sp_output_file.
	 const char*const name_part,      ///< Defined for \ref fopen_sp_output_file.
	 const char*const extension_part, ///< Defined for \ref fopen_sp_output_file.
	 const int mpi_rank               ///< Defined for \ref fopen_sp_output_file.
	);

/// \brief Set up the input file name for the given specifier if not already done.
static void set_up_input_name
	(const char input_spec,           ///< Defined for \ref fopen_input.
	 char*const input_name_spec,      ///< The specialized input name to be set up.
	 bool*const need_set_up,          ///< Flag for whether the name needs to be set up.
	 const char*const ctrl_name_full, ///< Defined for \ref fopen_input.
	 const char*const input_path      ///< Defined for \ref fopen_input.
	);

// Interface functions ********************************************************************************************** //

const char* extract_path (const char*const name_full)
{
	const size_t len = strlen(name_full);

	size_t path_len = len;
	while (name_full[path_len] != '/')
		--path_len; // Do nothing

	assert(path_len != len);
	assert(path_len != 0);

	static char name_path[STRLEN_MAX] = { 0, };
	for (size_t i = 0; i < path_len; ++i)
		name_path[i] = name_full[i];
	name_path[path_len] = 0;

	return name_path;
}

FILE* fopen_checked (const char*const file_name_full)
{
	FILE* file = fopen(file_name_full,"r");
	if (file == NULL) {
		printf("File: '%s' is not present.\n",file_name_full);
		EXIT_UNSUPPORTED;
	}
	return file;
}

void set_up_fopen_input (const char*const ctrl_name_full, const char*const input_path)
{
	FILE*const tmp = fopen_input(0,ctrl_name_full,input_path);
	assert(tmp == NULL);
}

FILE* fopen_input (const char input_spec, const char*const ctrl_name_full_i, const char*const input_path_i)
{
	static char ctrl_name_full[STRLEN_MAX] = {0};
	static char input_path[STRLEN_MAX]     = {0};
	static bool need_set_up = true;

	assert(input_spec == 0 || !need_set_up);

	char* input_name;
	switch (input_spec) {
	case 0: {
		// Only performs set up (Does not return a `FILE*`).
		need_set_up = false;
		strcpy(ctrl_name_full,ctrl_name_full_i);
		strcpy(input_path,input_path_i);
		return NULL;
		break;
	} case 'c': {
		input_name = ctrl_name_full;
		break;
	} case 'g': {
		static char input_name_spec[STRLEN_MAX] = {0};
		static bool need_set_up_name = true;
		set_up_input_name(input_spec,input_name_spec,&need_set_up_name,ctrl_name_full,input_path);

		input_name = input_name_spec;
		break;
	} case 's': {
		static char input_name_spec[STRLEN_MAX] = {0};
		static bool need_set_up_name = true;
		set_up_input_name(input_spec,input_name_spec,&need_set_up_name,ctrl_name_full,input_path);

		input_name = input_name_spec;
		break;
	} case 't': {
		static char input_name_spec[STRLEN_MAX] = {0};
		static bool need_set_up_name = true;
		set_up_input_name(input_spec,input_name_spec,&need_set_up_name,ctrl_name_full,input_path);

		// MSB: Return the input name found for the test case data file
		input_name = input_name_spec;
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",input_spec);
		break;
	}

	return fopen_checked(input_name);
}

FILE* fopen_create_dir (const char*const file_name_full)
{
	const char*const file_path = extract_path(file_name_full);
	mkdir_p(file_path);
	return fopen(file_name_full,"w");
}

FILE* fopen_sp_output_file
	(const char sp_type, const char*const name_part, const char*const extension_part, const int mpi_rank)
{
	const char* file_name = get_file_name_sp(sp_type,name_part,extension_part,mpi_rank);
	return fopen_create_dir(file_name);
}

FILE* fopen_sp_input_file
	(const char sp_type, const char*const name_part, const char*const extension_part, const int mpi_rank)
{
	const char* file_name = get_file_name_sp(sp_type,name_part,extension_part,mpi_rank);
	return fopen_checked(file_name);
}

FILE* fopen_sp_input_file_unchecked
	(const char sp_type, const char*const name_part, const char*const extension_part, const int mpi_rank)
{
	const char* file_name = get_file_name_sp(sp_type,name_part,extension_part,mpi_rank);
	return fopen(file_name,"r");
}

void mkdir_p_given_file_name (const char*const file_name)
{
	const char*const file_path = extract_path(file_name);
	mkdir_p(file_path);
}

void fgets_checked (char*const line, const int sizeof_line, FILE*const file)
{
	if (!fgets(line,sizeof_line,file))
		EXIT_ERROR("fgets returned NULL.\n");
}

bool first_string_matches (const char*const line, const char*const str_search)
{
	char str_first[STRLEN_MAX];
	sscanf(line,"%s",str_first);
	if (!strstr(line,str_search) || (strcmp(str_search,str_first) != 0))
		return false;
	return true;
}

void skip_lines (FILE* file, const int n_skip)
{
	char line[STRLEN_MAX];
	for (int n = 0; n < n_skip; ++n)
		fgets_checked(line,sizeof(line),file);
}

void skip_lines_ptr (FILE* file, char**const line, const int line_size, const int n_skip)
{
	for (int n = 0; n < n_skip; ++n)
		fgets_checked(*line,line_size,file);
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

void read_skip_name_const_c_1 (const char*const var_name, const char*const line, const char*const var)
{
	if (!first_string_matches(line,var_name))
		return;

	sscanf(line,"%*s %s",(char*)var);
}

void read_skip_c (const char*const line, char*const var)
{
	sscanf(line,"%*s %s",var);
}

void read_skip_name_i (const char*const var_name, const char*const line, int*const var)
{
	if (!first_string_matches(line,var_name))
		return;
	read_skip_i(line,var);
}

void read_skip_i (const char*const line, int*const var)
{
	sscanf(line,"%*s %u",var);
}

void read_skip_const_i (const char*const line, const int*const var)
{
	read_skip_i(line,(int*)var);
}

void read_skip_c_1 (const char*const line, char*const var)
{
	sscanf(line,"%*s %s",var);
}

void read_skip_const_c_1 (const char*const line, const char*const var)
{
	read_skip_c_1(line,(char*)var);
}

void read_skip_c_2 (const char*const line_i, const int n_skip, char**const var, const int n_var)
{
	char line[STRLEN_MAX];
	strcpy(line,line_i);

	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		sscanf(token_s,"%s",var[i]);
		token_s = strtok(NULL," ");
	}
}

void read_skip_const_b (const char*const line, const bool*const var)
{
	int tmp = 0;
	sscanf(line,"%*s %d",&tmp);
	const_cast_b(var,tmp);
}

void read_skip_d (char*const line, double*const var, const int n_skip, const bool remove_semi)
{
	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	char* token_s_clean = ( remove_semi ? strtok(token_s,";") : token_s );

	sscanf(token_s_clean,"%lf",var);
}

void read_skip_const_d (char*const line, const double*const var, const int n_skip, const bool remove_semi)
{
	read_skip_d(line,(double*)var,n_skip,remove_semi);
}

void read_skip_i_1 (char*const line_i, const int n_skip, int*const var, const int n_var)
{
	char line[STRLEN_MAX];
	strcpy(line,line_i);

	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		sscanf(token_s,"%d",&var[i]);
		token_s = strtok(NULL," ");
	}
}

void read_skip_const_i_1 (char*const line_i, const int n_skip, const int*const var, const int n_var)
{
	read_skip_i_1(line_i,n_skip,(int*)var,n_var);
}

void read_skip_ptrdiff_1 (char*const line_i, const int n_skip, ptrdiff_t*const var, const int n_var)
{
	char line[STRLEN_MAX];
	strcpy(line,line_i);

	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		sscanf(token_s,"%td",&var[i]);
		token_s = strtok(NULL," ");
	}
}

void read_skip_d_1 (char*const line_i, const int n_skip, double*const var, const int n_var)
{
	char line[STRLEN_MAX];
	strcpy(line,line_i);

	char* token_s = strtok(line," ");
	for (int i = 0; i < n_skip; ++i)
		token_s = strtok(NULL," ");

	for (int i = 0; i < n_var; ++i) {
		sscanf(token_s,"%lf",&var[i]);
		token_s = strtok(NULL," ");
	}
}

void read_skip_const_d_1 (char*const line_i, const int n_skip, const double*const var, const int n_var)
{
	read_skip_d_1(line_i,n_skip,(double*)var,n_var);
}

void read_skip_string_count_i (const char*const str_search, int*const count, char*const line_i, int*const var)
{
	if (!first_string_matches(line_i,str_search))
		return;

	(*count)++;
	read_skip_i(line_i,var);
}

void read_skip_string_count_d (const char*const str_search, int*const count, char*const line_i, double*const var)
{
	if (!first_string_matches(line_i,str_search))
		return;

	(*count)++;
	read_skip_d(line_i,var,1,false);
}

void read_skip_string_count_const_d
	(const char*const str_search, int*const count, char*const line_i, const double*const var)
{
	read_skip_string_count_d(str_search,count,line_i,(double*)var);
}

void read_skip_string_count_c_style_d
	(const char*const str_search, int*const count, char*const line_i, double*const var)
{
	if (!first_string_matches(line_i,str_search))
		return;

	(*count)++;
	read_skip_d(line_i,var,2,true);
}

void read_skip_string_count_c_style_const_d
	(const char*const str_search, int*const count, char*const line_i, const double*const var)
{
	read_skip_string_count_c_style_d(str_search,count,line_i,(double*)var);
}

void read_skip_file_const_b (const char*const var_name, FILE* file, const bool*const var)
{
	char line[STRLEN_MAX];
	fgets_checked(line,sizeof(line),file);

	if (!strstr(line,var_name))
		EXIT_ERROR("Did not find '%s' in the current line of the file.\n",var_name);

	read_skip_const_b(line,var);
}

void read_skip_file_const_i (const char*const var_name, FILE* file, const int*const var)
{
	char line[STRLEN_MAX];
	fgets_checked(line,sizeof(line),file);

	if (!strstr(line,var_name))
		EXIT_ERROR("Did not find '%s' in the current line of the file.\n",var_name);

	read_skip_const_i_1(line,1,var,1);
}

void read_skip_file_i (const char*const var_name, FILE* file, int*const var)
{
	char line[STRLEN_MAX];
	fgets_checked(line,sizeof(line),file);

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
	const char* beg = &name_full[0], // Pointer to the first element of the name_full.
	          * end = NULL;          // Pointer to one past the last element of the name_full.

	bool found_extension = false;
	int ind = 0;
	for (int i = (int)strlen(name_full)-1; i >= 0; --i) {
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

	static char name[STRLEN_MAX] = { 0, };

	ind = 0;
	while (end != beg) {
		name[ind] = *beg++;
		++ind;
	}
	name[ind] = 0;

	return name;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Converts the string pointed to by the argument `str` to a smaller type than `long int`.
 *  \return See brief.
 *
 *  The implementation was taken from [this SO answer][SO_strtoi].
 *
 *  <!-- References: -->
 *  [SO_strtoi]: https://stackoverflow.com/a/29378380/5983549
 */
static long strto_subrange
	(const char *str, ///< See strtoi.
	 char **endptr,   ///< See strtoi.
	 int base,        ///< See strtoi.
	 long min,        ///< The minimum supported value.
	 long max         ///< The maximum supported value.
	);

/** \brief Search for an extension in the control file with the given input name.
 *  \return A static `char*` holding the name of the extension if found, `NULL` otherwise. */
static const char* read_extension
	(const char*const ctrl_name_full, ///< The filename of the file in which to search for the extension.
	 const char*const extension_name  ///< The name of the extension to search for.
	);

static int strtoi (const char* str, char** endptr, int base)
{
#if INT_MAX == LONG_MAX && INT_MIN == LONG_MIN
	return (int) strtol(str,endptr,base);
#else
	return (int) strto_subrange(str,endptr,base,INT_MIN,INT_MAX);
#endif
}

static int mkdir_p (const char *path)
{
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */
    const size_t len = strlen(path);
    char _path[PATH_MAX];
    char *p;

    errno = 0;

    /* Copy string so its mutable */
    if (len > sizeof(_path)-1) {
        errno = ENAMETOOLONG;
        return -1;
    }
    strcpy(_path, path);

    /* Iterate the string */
    for (p = _path + 1; *p; p++) {
        if (*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if (mkdir(_path, S_IRWXU) != 0) {
                if (errno != EEXIST)
                    return -1;
            }

            *p = '/';
        }
    }

    if (mkdir(_path, S_IRWXU) != 0) {
        if (errno != EEXIST)
            return -1;
    }

    return 0;
}

static const char* get_file_name_sp
	(const char sp_type, const char*const name_part, const char*const extension_part, const int mpi_rank)
{
	assert(sp_type == 's' || sp_type == 'p');
	static char file_name[STRLEN_MAX] = { 0, };

	if (sp_type == 's') {
		sprintf(file_name,"%s%c%d%s%s",name_part,'_',mpi_rank,".",extension_part);
	} else if (sp_type == 'p') {
		if (strcmp(extension_part,"txt") == 0)
			sprintf(file_name,"%s%s%s",name_part,"_p.",extension_part);
		else
			sprintf(file_name,"%s%s%s",name_part,".p",extension_part);
	}
	return file_name;
}

static void set_up_input_name
	(const char input_spec, char*const input_name_spec, bool*const need_set_up, const char*const ctrl_name_full,
	 const char*const input_path)
{
	if (!*need_set_up)
		return;
	*need_set_up = false;

	int index = sprintf(input_name_spec,"%s",input_path);
	switch (input_spec) {
	case 'g': {
		const char*const extension = read_extension(ctrl_name_full,"geometry_parameters_extension");
		if (extension[0] == 0)
			sprintf(index+input_name_spec,"%s","geometry_parameters.geo");
		else
			sprintf(index+input_name_spec,"%s%s%s","geometry_parameters_",extension,".geo");
		break;
	} case 's': {
		const char*const extension = read_extension(ctrl_name_full,"solution_extension");
		if (extension[0] == 0)
			sprintf(index+input_name_spec,"%s","solution.data");
		else
			sprintf(index+input_name_spec,"%s%s%s","solution_",extension,".data");
		break;
	} case 't': {
		// MSB: Set the test case absolute file path and return it
		const char*const extension = read_extension(ctrl_name_full,"test_case_extension");
		// MSB: In our test case, we have no test_case_extension parameter so we will just read
		// the test_case.data value (extension[0] will == 0)
		if (extension[0] == 0)
			sprintf(index+input_name_spec,"%s","test_case.data");
		else
			sprintf(index+input_name_spec,"%s%s%s","test_case_",extension,".data");
		break;
	} default:
		EXIT_ERROR("Unsupported: %c\n",input_spec);
		break;
	}
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

static const char* read_extension (const char*const ctrl_name_full, const char*const extension_name)
{
	static char extension[STRLEN_MAX];
	extension[0] = 0;

	FILE*const ctrl_file = fopen_checked(ctrl_name_full);
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),ctrl_file))
		read_skip_name_const_c_1(extension_name,line,extension);
	fclose(ctrl_file);

	return extension;
}
