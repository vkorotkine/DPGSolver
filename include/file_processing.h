// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__file_processing_h__INCLUDED
#define DPG__file_processing_h__INCLUDED
/**	\file
 *	Provides file processing related functions.
 */

#include <stdio.h>
#include <stdbool.h>

// Opening files **************************************************************************************************** //

/// \brief Open file and check for successful completion.
FILE* fopen_checked
	(const char*const file_name_full ///< File name including full path.
	);

// Reading data from the current line ******************************************************************************* //

/// \brief Skip lines in a file while reading
void skip_lines
	(FILE* file,               ///< Standard.
	 char**const line,         ///< The line.
	 const size_t line_size,   ///< The size allocated for the line.
	 const unsigned int n_skip ///< The number of lines to skip.
	);

/// \brief Discard values from the beginning of a line.
void discard_line_values
	(char**const line,      ///< The line.
	 unsigned int n_discard ///< The number of values to discard.
	);

/// \brief Reads values from the line into an `unsigned int` array.
void read_line_values_ui
	(char**const line,         ///< The line.
	 const unsigned int n_val, ///< The number of values to read.
	 unsigned int*const vals,  ///< The array in which to store the values.
	 const bool decrement      /**< Flag for whether decrementing by 1 is enabled. Used to convert from 1-based to
	                            *   0-based indexing. */
	);

/// \brief Read a `char*`, skipping the first string.
void read_skip_c
	(const char*const line, ///< Line from which to read data.
	 char*const var         ///< Variable in which to store data.
	);

/// \brief Read an `unsigned int*`, skipping the first string.
void read_skip_ui
	(const char*const line, ///< Line from which to read data.
	 unsigned int*const var ///< Variable in which to store data.
	);

/// \brief Read a `const char*`, skipping the first string.
void read_skip_const_c
	(const char*const line, ///< Line from which to read data.
	 const char*const var   ///< Variable in which to store data.
	);

/// \brief Read a `const unsigned int*`, skipping the first string.
void read_skip_const_ui
	(const char*const line,        ///< Line from which to read data.
	 const unsigned int *const var ///< Variable in which to store data.
	);

/// \brief Read a `const bool*`, skipping the first string.
void read_skip_const_b
	(const char*const line, ///< Line from which to read data.
	 const bool *const var  ///< Variable in which to store data.
	);

// Setting file names with paths ************************************************************************************ //

/// \brief Append `src` (char*) to `dest` with optional forward slash ('\').
void strcat_path_c
	(char* dest,           ///< Destination.
	 const char*const src, ///< Source.
	 bool add_slash        ///< Flag indicating whether a forward slash ('\') should be appended.
	);

/// \brief Append `src` (unsigned int) to `dest`.
void strcat_path_ui
	(char* dest,            ///< Destination.
	 const unsigned int src ///< Source.
	);

#endif // DPG__file_processing_h__INCLUDED
