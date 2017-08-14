// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Element_h__INCLUDED
#define DPG__Element_h__INCLUDED
/**	\file
 *
 *	\todo Add comments.
 */

#include "Intrusive.h"

/// \brief Container for data relating to the base Elements.
struct Element {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const unsigned int d; ///< The dimension.

	/// The (corr)espondence between the (f)aces and (ve)rtices.
	const struct const_Multiarray_Vector_ui*const corr_f_ve;

};

/// \brief Constructs the base \ref Element \ref Intrusive_List.
struct Intrusive_List* constructor_Element_List
	(const unsigned int d ///< The dimension.
	);

#endif // DPG__Element_h__INCLUDED
