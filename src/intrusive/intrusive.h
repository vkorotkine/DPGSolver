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
 *
 *  While not required for the conventional usage of the intrusive list functionality, certain members have been added
 *  to the \ref Intrusive_List and \ref Intrusive_Link containers to facilitate the use of lists of derived objects from
 *  lists of base objects:
 *  - \ref Intrusive_List::base provides a pointer to the base list (if applicable) such that the base list can be
 *    destructed after all pointers to its links have been changed to pointers to derived links.
 *  - \ref Intrusive_Link::derived provides a pointer to the derived link (if applicable) such that all pointers to base
 *    links can be replaced with pointers to derived links.
 */

#include <stdbool.h>

/// \brief A doubly-linked list structure to hold intrusive containers.
struct Intrusive_List {
	struct Intrusive_Link* first; ///< Pointer to the first \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* last;  ///< Pointer to the last \ref Intrusive_Link\* in the list.

	int name; ///< The name of the list. Used for selecting the appropriate destructor function.
	struct Intrusive_List* base;  ///< Pointer to the base \ref Intrusive_List\* if applicable.
};

/// \brief `const` version of \ref Intrusive_List.
struct const_Intrusive_List {
	const struct const_Intrusive_Link*const first; ///< Defined in \ref Intrusive_List.
	const struct const_Intrusive_Link*const last;  ///< Defined in \ref Intrusive_List.

	const int name; ///< Defined in \ref Intrusive_List.
	const struct const_Intrusive_List*const base;  ///< Defined in \ref Intrusive_List.
};

/// \brief A link for a doubly-linked list.
struct Intrusive_Link {
	struct Intrusive_Link* prev; ///< Pointer to the previous \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* next; ///< Pointer to the next \ref Intrusive_Link\* in the list.

	struct Intrusive_Link* derived; ///< Pointer to the derived \ref Intrusive_Link\* if applicable.
};

/** \brief `const` version of \ref Intrusive_Link.
 *  \note Pointers are not `const` as this precludes iterating over them.
 */
struct const_Intrusive_Link {
	const struct const_Intrusive_Link* prev; ///< Defined in \ref Intrusive_Link.
	const struct const_Intrusive_Link* next; ///< Defined in \ref Intrusive_Link.

	const struct const_Intrusive_Link*const derived; ///< Defined in \ref Intrusive_Link.
};

// Interface functions ********************************************************************************************** //

/** \brief Contructs an empty \ref Intrusive_List.
 *  \return Standard. */
struct Intrusive_List* constructor_empty_IL
	(const int list_name,        ///< \ref Intrusive_List::name.
	 struct Intrusive_List* base ///< \ref Intrusive_List::base.
	);

/** \brief `const` version of \ref constructor_empty_IL.
 *  \return Standard. */
const struct const_Intrusive_List* constructor_empty_const_IL
	(const int list_name,                    ///< Defined for \ref constructor_empty_IL.
	 const struct const_Intrusive_List* base ///< Defined for \ref constructor_empty_IL.
	);

/** \brief Destructs all Intrusive_Links in the \ref Intrusive_List.
 *  \note Only frees the link entries and not any members which may be pointing to allocated memory. */
void clear_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/** \brief Destructs a \ref Intrusive_List and optionally its dynamically allocated \ref Intrusive_Link components.
 *
 *  \note No destructor is being called for the elements represented by the \ref Intrusive_Link components.
 */
void destructor_IL
	(struct Intrusive_List* lst, ///< Standard.
	 const bool destruct_links   ///< Flag for whether the individual links require destruction.
	);

/// \brief `const` version of \ref destructor_IL.
void destructor_const_IL
	(const struct const_Intrusive_List* lst, ///< See brief.
	 const bool destruct_links               ///< See brief.
	);

/// \brief Destructor for a \ref Intrusive_List which is a base for the input \ref Intrusive_List.
void destructor_IL_base
	(struct Intrusive_List* lst ///< Standard.
	);

/// \brief `const` version of \ref destructor_IL_base.
void destructor_const_IL_base
	(const struct const_Intrusive_List* lst ///< Defined for \ref destructor_IL_base.
	);

/// \brief Add an \ref Intrusive_Link to the end of the \ref Intrusive_List.
void push_back_IL
	(struct Intrusive_List* lst, ///< The list.
	 struct Intrusive_Link* curr ///< The current link.
	);

/// \brief `const` version of \ref push_back_IL.
void push_back_const_IL
	(const struct const_Intrusive_List* lst, ///< Defined for \ref push_back_IL.
	 const struct const_Intrusive_Link* curr ///< Defined for \ref push_back_IL.
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

/** \brief `const` version of \ref erase_IL.
 *  \return See brief. */
const struct const_Intrusive_Link* erase_const_IL
	(const struct const_Intrusive_List* lst,  ///< Defined for \ref erase_IL.
	 const struct const_Intrusive_Link* curr  ///< Defined for \ref erase_IL.
	);

/** \brief Replace the input number of links in the list with those of the new list.
 *
 *  \warning Does not free memory.
 *
 *  \return Pointer to the last link of the new list elements.
 */
struct Intrusive_Link* replace_IL
	(struct Intrusive_List*const lst,    ///< The list.
	 const int n_replace,                ///< The number of links in the list to be erased.
	 struct Intrusive_Link*const curr_f, ///< The current link (first in group to be replaced).
	 struct Intrusive_Link*const first   ///< The first link of the new sub-list.
	);

/** \brief Return the \ref Intrusive_Link `n_adv` links down from the input link.
 *  \return See brief. */
struct Intrusive_Link* advance_Link
	(const int n_adv,             ///< The number of links to advance.
	 struct Intrusive_Link* first ///< The first link.
	);

/// \brief Insert a "sub" list into the "main" list before the "curr"ent link.
void insert_List_into_List
	(struct Intrusive_List*const sub,  ///< The sub list.
	 struct Intrusive_List*const main, ///< The main list.
	 struct Intrusive_Link*const curr  ///< The current link.
	);

#endif // DPG__Intrusive_h__INCLUDED
