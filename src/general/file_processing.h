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

#ifndef DPG__file_processing_h__INCLUDED
#define DPG__file_processing_h__INCLUDED
/** \file
 *  \brief Provides file processing related functions.
 */

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

#include "file_processing_conversions.h"

// Opening files **************************************************************************************************** //

/** \brief Extract the path from the string input.
 *  \return The path (no free necessary). */
const char* extract_path
	(const char*const name_full ///< The full input.
	);

/** \brief Open file and check for successful completion.
 *  \return See brief. */
FILE* fopen_checked
	(const char*const file_name_full ///< File name including full path.
	);

/// \brief Set up the static variables in \ref fopen_input.
void set_up_fopen_input
	(const char*const ctrl_name_full, ///< Full path to a control file in which to filename extensions may be found.
	 const char*const input_path      ///< Full path to input file directory.
	);

/** \brief Open the an input file based on the input parameters.
 *  \return See brief. */
FILE* fopen_input
	(const char input_spec,             ///< The input specifier. Options: 'c'ontrol, 'g'eometry, 's'olution, 't'est case, 'o'ptimization.
	 const char*const ctrl_name_full_i, ///< Defined for \ref set_up_fopen_input.
	 const char*const input_path_i      ///< Defined for \ref set_up_fopen_input.
	);

/** \brief Open file of the given input name, creating the directory if it does not exist.
 *  \return See brief. */
FILE* fopen_create_dir
	(const char*const file_name_full ///< File name including full path.
	);

/** \brief Open the file to which the 's'erial/'p'arallel output will be written.
 *  \return The pointer to the file. */
FILE* fopen_sp_output_file
	(const char sp_type,              ///< Type indicator for 's'erial or 'p'arallel.
	 const char*const name_part,      ///< Partial output file name.
	 const char*const extension_part, ///< File extension.
	 const int mpi_rank               ///< The mpi rank.
	);

/** \brief Open the file from which the 's'erial/'p'arallel input will be read.
 *  \return The pointer to the file. */
FILE* fopen_sp_input_file
	(const char sp_type,              ///< Defined for \ref fopen_sp_output_file.
	 const char*const name_part,      ///< Defined for \ref fopen_sp_output_file.
	 const char*const extension_part, ///< Defined for \ref fopen_sp_output_file.
	 const int mpi_rank               ///< Defined for \ref fopen_sp_output_file.
	);

/** \brief Unchecked version of \ref fopen_sp_input_file.
 *  \return The pointer to the file or `NULL` if not present. */
FILE* fopen_sp_input_file_unchecked
	(const char sp_type,              ///< See brief.
	 const char*const name_part,      ///< See brief.
	 const char*const extension_part, ///< See brief.
	 const int mpi_rank               ///< See brief.
	);

/// \brief Create a directory associated with the path of the input name if not already present.
void mkdir_p_given_file_name
	(const char*const file_name ///< The name of the file containing the path to the directory.
	);

// Reading data from the current line ******************************************************************************* //

/// \brief Call standard library fgets, exiting if `NULL` is returned.
void fgets_checked
	(char*const line,       ///< Input[0] for fgets.
	 const int sizeof_line, ///< Input[1] for fgets.
	 FILE*const file        ///< Input[2] for fgets.
	);

/** \brief Check whether the first string in the line is equal to the desired input string.
 *  \return `true` if yes; `false` otherwise. */
bool first_string_matches
	(const char*const line,      ///< The line to be checked.
	 const char*const str_search ///< The string to search for as the first entry of the line.
	);

/// \brief Skip lines in a file while reading.
void skip_lines
	(FILE* file,          ///< Standard.
	 const int n_skip     ///< The number of lines to skip.
	);

/// \brief Skip lines in a file while reading, supporting a pointer to the next line.
void skip_lines_ptr
	(FILE* file,          ///< Standard.
	 char**const line,    ///< The line.
	 const int line_size, ///< The size allocated for the line.
	 const int n_skip     ///< The number of lines to skip.
	);

/// \brief Discard values from the beginning of a line.
void discard_line_values
	(char**const line, ///< The line.
	 int n_discard     ///< The number of values to discard.
	);

/// \brief Reads values from the line into an `int` array.
void read_line_values_i
	(char**const line,      ///< The line.
	 const ptrdiff_t n_val, ///< The number of values to read.
	 int*const vals,        ///< The array in which to store the values.
	 const bool decrement   /**< Flag for whether decrementing by 1 is enabled. Used to convert from 1-based to
	                         *   0-based indexing. */
	);

/** \brief Reads values from the line into a `long int` array.
 *  \note This function is also currently used to read values into `ptrdiff_t` arrays as there is not standard library
 *        function associated with this data type. */
void read_line_values_l
	(char**const line,      ///< The line.
	 const ptrdiff_t n_val, ///< The number of values to read.
	 long int*const vals,   ///< The array in which to store the values.
	 const bool decrement   /**< Flag for whether decrementing by 1 is enabled. Used to convert from 1-based to
	                         *   0-based indexing. */
	);

/// \brief Reads values from the line into a `double` array.
void read_line_values_d
	(char**const line,      ///< The line.
	 const ptrdiff_t n_val, ///< The number of values to read.
	 double*const vals      ///< The array in which to store the values.
	);

/** \brief Read the second string in the line into a `const char*` if the first string in the line matches the input
 *         variable name. */
void read_skip_name_const_c_1
	(const char*const var_name, ///< The variable name for which to search.
	 const char*const line,     ///< Line from which to read data.
	 const char*const var       ///< Variable in which to store data.
	);

/// \brief Read a `char`, skipping the first string.
void read_skip_c
	(const char*const line, ///< Line from which to read data.
	 char*const var         ///< Variable in which to store data.
	);

/** \brief Read the second string in the line into a `int*` if the first string in the line matches the input
 *         variable name. */
void read_skip_name_i
	(const char*const var_name, ///< The variable name for which to search.
	 const char*const line,     ///< Line from which to read data.
	 int*const var              ///< Variable in which to store data.
	);

/// \brief Read an `int`, skipping the first string.
void read_skip_i
	(const char*const line, ///< Line from which to read data.
	 int*const var          ///< Variable in which to store data.
	);

/// \brief `const` version of \ref read_skip_i.
void read_skip_const_i
	(const char*const line, ///< Defined for \ref read_skip_i.
	 const int*const var    ///< Defined for \ref read_skip_i.
	);

/// \brief Read a `char*`, skipping the first string.
void read_skip_c_1
	(const char*const line, ///< Line from which to read data.
	 char*const var         ///< Variable in which to store data.
	);

/// \brief `const` version of \ref read_skip_c_1.
void read_skip_const_c_1
	(const char*const line, ///< See brief.
	 const char*const var   ///< See brief.
	);

/// \brief Read into a preallocated `char**`, optionally skipping strings.
void read_skip_c_2
	(const char*const line_i, ///< Line from which to read data.
	 const int n_skip,        ///< The number of strings to skip.
	 char**const var,         ///< Variable in which to store data.
	 const int n_var          ///< The number of entries to store in the `var` array.
	);

/// \brief Read a `bool`, skipping the first string.
void read_skip_b
	(const char*const line, ///< Line from which to read data.
	 bool*const var         ///< Variable in which to store data.
	);

/// \brief `const` version of \ref read_skip_b.
void read_skip_const_b
	(const char*const line, ///< See brief.
	 const bool*const var   ///< See brief.
	);

/// \brief Read a `double`, optionally skipping strings and optionally removing trailing semicolons.
void read_skip_d
	(char*const line,       ///< Line from which to read data.
	 double*const var,      ///< Variable in which to store data.
	 const int n_skip,      ///< The number of strings to skip.
	 const bool remove_semi ///< Flag for optional removal of semicolon.
	);

/// \brief `const` version of \ref read_skip_d.
void read_skip_const_d
	(char*const line,        ///< Defined for \ref read_skip_d.
	 const double*const var, ///< Defined for \ref read_skip_d.
	 const int n_skip,       ///< Defined for \ref read_skip_d.
	 const bool remove_semi  ///< Defined for \ref read_skip_d.
	);

/// \brief Get the next line from the input file and read the `var_name` variable into `var` of type `const bool`.
void read_skip_file_const_b
	(const char*const var_name, ///< The name of the Variable to search for.
	 FILE* file,                ///< File from which to read data.
	 const bool*const var       ///< Variable in which to store data.
	);

/// \brief Get the next line from the input file and read the `var_name` variable into `var` of type `const int`.
void read_skip_file_const_i
	(const char*const var_name, ///< The name of the Variable to search for.
	 FILE* file,                ///< File from which to read data.
	 const int*const var        ///< Variable in which to store data.
	);

/// \brief Get the next line from the input file and read the `var_name` variable into `var` of type `int`.
void read_skip_file_i
	(const char*const var_name, ///< The name of the Variable to search for.
	 FILE* file,                ///< File from which to read data.
	 int*const var              ///< Variable in which to store data.
	);

/// \brief Read a `int*`, optionally skipping strings.
void read_skip_i_1
	(char*const line_i, ///< Line from which to read data.
	 const int n_skip,  ///< The number of strings to skip.
	 int*const var,     ///< Variable in which to store data.
	 const int n_var    ///< The number of entries to store in the `var` array.
	);

/// \brief `const` version of \ref read_skip_i_1.
void read_skip_const_i_1
	(char*const line_i,   ///< Defined for \ref read_skip_i_1.
	 const int n_skip,    ///< Defined for \ref read_skip_i_1.
	 const int*const var, ///< Defined for \ref read_skip_i_1.
	 const int n_var      ///< Defined for \ref read_skip_i_1.
	);

/// \brief Read a `ptrdiff_t*`, optionally skipping strings.
void read_skip_ptrdiff_1
	(char*const line_i,   ///< Line from which to read data.
	 const int n_skip,    ///< The number of strings to skip.
	 ptrdiff_t*const var, ///< Variable in which to store data.
	 const int n_var      ///< The number of entries to store in the `var` array.
	);

/// \brief Read a `double*`, optionally skipping strings.
void read_skip_d_1
	(char*const line_i, ///< Line from which to read data.
	 const int n_skip,  ///< The number of strings to skip.
	 double*const var,  ///< Variable in which to store data.
	 const int n_var    ///< The number of entries to store in the `var` array.
	);

/// \brief `const` version of \ref read_skip_d_1.
void read_skip_const_d_1
	(char*const line_i,      ///< Defined for \ref read_skip_d_1.
	 const int n_skip,       ///< Defined for \ref read_skip_d_1.
	 const double*const var, ///< Defined for \ref read_skip_d_1.
	 const int n_var         ///< Defined for \ref read_skip_d_1.
	);

/// \brief Read an `int` from the input line and increment the count if the search string is in the line.
void read_skip_string_count_i
	(const char*const str_search, ///< String to search for in the line.
	 int*const count,             ///< Count to increment.
	 char*const line_i,           ///< Line from which to read data.
	 int*const var                ///< Variable in which to store the data.
	);

/// \brief Read a `double` from the input line and increment the count if the search string is in the line.
void read_skip_string_count_d
	(const char*const str_search, ///< String to search for in the line.
	 int*const count,             ///< Count to increment.
	 char*const line_i,           ///< Line from which to read data.
	 double*const var             ///< Variable in which to store the data.
	);

/// \brief `const` version of \ref read_skip_string_count_d.
void read_skip_string_count_const_d
	(const char*const str_search, ///< See brief.
	 int*const count,             ///< See brief.
	 char*const line_i,           ///< See brief.
	 const double*const var       ///< See brief.
	);

/** \brief Read a `double` from the input line and increment the count if the search string is in the line.
 *
 *  This function assumes that the input line has the form of a c-style variable definition (i.e. var = val;) and thus
 *  skips the first 2 entries in the string and removes the trailing semi-colon.
 */
void read_skip_string_count_c_style_d
	(const char*const str_search, ///< String to search for in the line.
	 int*const count,             ///< Count to increment.
	 char*const line_i,           ///< Line from which to read data.
	 double*const var             ///< Variable in which to store the data.
	);

/// \brief `const` version of \ref read_skip_string_count_c_style_d.
void read_skip_string_count_c_style_const_d
	(const char*const str_search, ///< See brief.
	 int*const count,             ///< See brief.
	 char*const line_i,           ///< See brief.
	 const double*const var       ///< See brief.
	);

// Setting/Getting file names *************************************************************************************** //

/// \brief Append `src` (char*) to `dest` with optional forward slash ('\').
void strcat_path_c
	(char* dest,            ///< Destination.
	 const char*const src,  ///< Source.
	 const char*const trail ///< Optional trailing characters to append.
	);

/// \brief Append `src` (int) to `dest`.
void strcat_path_i
	(char* dest,   ///< Destination.
	 const int src ///< Source.
	);

/** \brief Extract the name of the file from the input string (i.e. excluding the path and the extension).
 *  \return See brief. */
char* extract_name
	(const char*const name_full,  ///< The full name.
	 const bool extension_present ///< Flag for whether the extension is present.
	);

#endif // DPG__file_processing_h__INCLUDED
