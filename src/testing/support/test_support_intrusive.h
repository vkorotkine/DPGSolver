// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_support_intrusive_h__INCLUDED
#define DPG__test_support_intrusive_h__INCLUDED
/**	\file
 *	\brief Provides support functions for testing relating to containers defined in \ref intrusive.h.
 */

struct Intrusive_List;
struct const_Intrusive_List;

/**	\brief Constructor for an intrusive list of type specified by the `list_name`.
 *	\return Standard. */
struct Intrusive_List* constructor_file_name_IL
	(const char*const list_name,                       ///< The name of the type of container to form the links.
	 const char*const file_name,                       ///< The name of the file from which to read the data.
	 const struct const_Intrusive_List*const elements, ///< \ref Simulation::elements.
	 const struct Intrusive_List*const volumes         ///< \ref Simulation::volumes.
	);

#endif // DPG__test_support_intrusive_h__INCLUDED
