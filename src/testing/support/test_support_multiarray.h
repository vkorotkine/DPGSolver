// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_multiarray_h__INCLUDED
#define DPG__test_support_multiarray_h__INCLUDED
/**	\file
 *	Provides support functions for testing relating to the containers defined in \ref multiarray.h.
 */

/** \brief Constructor for a \ref Multiarray_Vector_i\* from data in the input file.
 *	\return Standard. */
struct Multiarray_Vector_i* constructor_file_Multiarray_Vector_i
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

#endif // DPG__test_support_multiarray_h__INCLUDED
