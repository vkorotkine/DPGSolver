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

#ifndef DPG__Intrusive_h__INCLUDED
#define DPG__Intrusive_h__INCLUDED
/** \file
 *  \brief Provides structures/functions related to the intrusive list functionality.
 *
 *  The data structures used here were inspired by the intrusive container example from ch. 27.9 of Stroustrup
 *  \cite Stroustrup2014. The motivating principle is that elements of an \ref Intrusive_List may be manipulated without
 *  knowing anything about the internal structure of the \ref Intrusive_Link.
 */

/// \brief A doubly-linked list structure to hold intrusive containers.
struct Intrusive_List {
	struct Intrusive_Link* first; ///< Pointer to the first \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* last;  ///< Pointer to the last \ref Intrusive_Link\* in the list.

	int name; ///< The name of the list. Used for selecting the appropriate destructor function.
};

/// \brief `const` version of \ref Intrusive_List.
struct const_Intrusive_List {
	const struct Intrusive_Link*const first; ///< Defined in \ref Intrusive_List.
	const struct Intrusive_Link*const last;  ///< Defined in \ref Intrusive_List.

	int name; ///< Defined in \ref Intrusive_List.
};

/// \brief A link for a doubly-linked list.
struct Intrusive_Link {
	struct Intrusive_Link* prev; ///< Pointer to the previous \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* next; ///< Pointer to the next \ref Intrusive_Link\* in the list.
};

/// \brief `const` version of \ref Intrusive_Link.
struct const_Intrusive_Link {
	const struct Intrusive_Link*const prev; ///< See non-`const`.
	const struct Intrusive_Link*const next; ///< See non-`const`.
};

/** \brief Contructs an empty \ref Intrusive_List.
 *	\return Standard. */
struct Intrusive_List* constructor_empty_IL
	(const int list_name ///< \ref Intrusive_List::name.
	);

/** \brief Destructs all Intrusive_Links in the \ref Intrusive_List.
 *	\note Only frees the link entries and not any members which may be pointing to allocated memory. */
void clear_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/** \brief Destructs a \ref Intrusive_List and its dynamically allocated \ref Intrusive_Link components.
 *
 *  \note No destructor is being called for the elements represented by the \ref Intrusive_Link components.
 */
void destructor_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/// \brief `const` version of \ref destructor_IL.
void destructor_const_IL
	(const struct const_Intrusive_List* lst ///< Standard.
	);

/// \brief Add an \ref Intrusive_Link to the end of the \ref Intrusive_List.
void push_back_IL
	(struct Intrusive_List* lst, ///< The list.
	 struct Intrusive_Link* curr ///< The current link.
	);

/** \brief Erase the current \ref Intrusive_Link.
 *
 *  \warning Does not free memory.
 *
 *  \return Pointer to the next \ref Intrusive_Link.
 */
struct Intrusive_Link* erase_IL
	(struct Intrusive_List* lst, ///< The list.
	 struct Intrusive_Link* curr ///< The current link.
	);

#endif // DPG__Intrusive_h__INCLUDED
