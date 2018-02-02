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

#ifndef DPG__file_processing_conversions_h__INCLUDED
#define DPG__file_processing_conversions_h__INCLUDED
/** \file
 *  \brief Provides file processing related functions used for conversion of input string inputs to macro definitions
 *         used throughout the code.
 */

/** \brief Read a string skipping the first string and set the value of the corresponding `const int` according to the
 *         available defitions for the input definition type. */
void read_skip_convert_i
	(const char*const line,     ///< Line from which to read data.
	 const char*const def_type, ///< The type of definition from which to check for a suitable value for `var`.
	 int*const var,             ///< Variable in which to store data.
	 int*const count_found      ///< Counter to be incremented if the variable is found. Pass `NULL` if not needed.
	);

/// \brief `const` version of \ref read_skip_convert_i.
void read_skip_convert_const_i
	(const char*const line,     ///< See brief.
	 const char*const def_type, ///< See brief.
	 const int*const var,       ///< See brief.
	 int*const count_found      ///< See brief.
	);

#endif // DPG__file_processing_conversions_h__INCLUDED
