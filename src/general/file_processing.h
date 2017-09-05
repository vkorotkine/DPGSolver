// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__file_processing_h__INCLUDED
#define DPG__file_processing_h__INCLUDED
/**	\file
 *	\brief Provides file processing related functions.
 */

#include <stdio.h>
#include <stdbool.h>
#include <stddef.h>

// Opening files **************************************************************************************************** //

/** \brief Open file and check for successful completion.
 *	\return See brief. */
FILE* fopen_checked
	(const char*const file_name_full ///< File name including full path.
	);

/** \brief Open the an input file based on the input parameters.
 *	\return See brief. */
FILE* fopen_input
	(const char*const input_path, ///< Full path to the location of the input file.
	 const char*const input_spec  ///< The input specifier. Options: "geometry".
	);

// Reading data from the current line ******************************************************************************* //

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
	 const bool decrement   /**< Flag for whether decrementing by 1 is enabled. Used to convert from 1-based to 0-based
	                         *   indexing. */
	);

/** \brief Reads values from the line into a `long int` array.
 *	\note This function is also currently used to read values into `ptrdiff_t` arrays as there is not standard library
 *	      function associated with this data type. */
void read_line_values_l
	(char**const line,      ///< The line.
	 const ptrdiff_t n_val, ///< The number of values to read.
	 long int*const vals,   ///< The array in which to store the values.
	 const bool decrement   /**< Flag for whether decrementing by 1 is enabled. Used to convert from 1-based to 0-based
	                         *   indexing. */
	);

/// \brief Reads values from the line into a `double` array.
void read_line_values_d
	(char**const line,      ///< The line.
	 const ptrdiff_t n_val, ///< The number of values to read.
	 double*const vals      ///< The array in which to store the values.
	);

/// \brief Read a `char`, skipping the first string.
void read_skip_c
	(const char*const line, ///< Line from which to read data.
	 char*const var         ///< Variable in which to store data.
	);

/// \brief Read an `int`, skipping the first string.
void read_skip_i
	(const char*const line, ///< Line from which to read data.
	 int*const var          ///< Variable in which to store data.
	);

/// \brief Read a `const char`, skipping the first string.
void read_skip_const_c
	(const char*const line, ///< Line from which to read data.
	 const char*const var   ///< Variable in which to store data.
	);

/// \brief Read a `const int`, skipping the first string.
void read_skip_const_i
	(const char*const line, ///< Line from which to read data.
	 const int*const var    ///< Variable in which to store data.
	);

/// \brief Read a `const bool`, skipping the first string.
void read_skip_const_b
	(const char*const line, ///< Line from which to read data.
	 const bool*const var   ///< Variable in which to store data.
	);

/// \brief Read a `const double`, optionally skipping strings and optionally removing trailing semicolons.
void read_skip_const_d
	(const char*const line,  ///< Line from which to read data.
	 const double*const var, ///< Variable in which to store data.
	 const int n_skip,       ///< The number of strings to skip.
	 const bool remove_semi  ///< Flag for optional removal of semicolon.
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

/// \brief Read a `ptrdiff_t*`, optionally skipping strings.
void read_skip_ptrdiff_1
	(char*const line,     ///< Line from which to read data.
	 const int n_skip,    ///< The number of strings to skip.
	 ptrdiff_t*const var, ///< Variable in which to store data.
	 const int n_var      ///< The number of entries to store in the `var` array.
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
 *	\return See brief. */
char* extract_name
	(const char*const name_full,  ///< The full name.
	 const bool extension_present ///< Flag for whether the extension is present.
	);

#endif // DPG__file_processing_h__INCLUDED
