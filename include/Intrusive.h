// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Intrusive_h__INCLUDED
#define DPG__Intrusive_h__INCLUDED
/**	\file
 *
 *	The data structures used here were inspired by the intrusive container example from ch. 27.9 of Stroustrup
 *	\cite Stroustrup2014. The motivating principle is that elements of an \ref Intrusive_List may be manipulated without
 *	knowing anything about the internal structure of the \ref Intrusive_Link.
 *
 *	\todo Add comments.
 */

/// \brief A doubly-linked list structure to hold intrusive containers.
struct Intrusive_List {
	struct Intrusive_Link* first; ///< Pointer to the first \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* last;  ///< Pointer to the last \ref Intrusive_Link\* in the list.
};

/// \brief A link for a doubly-linked list.
struct Intrusive_Link {
	struct Intrusive_Link* prev; ///< Pointer to the previous \ref Intrusive_Link\* in the list.
	struct Intrusive_Link* next; ///< Pointer to the next \ref Intrusive_Link\* in the list.
};

/// \brief Initialize an \ref Intrusive_List to `NULL`.
void init_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/** \brief Contructs an empty \ref Intrusive_List.
 *	\return Pointer to the \ref Intrusive_List.
 */
struct Intrusive_List* constructor_empty_IL ();

/// \brief Destructs all Intrusive_Links in the \ref Intrusive_List.
void clear_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/** \brief Destructs a \ref Intrusive_List and its dynamically allocated \ref Intrusive_Link components.
 *
 *	\note No destructor is being called for the elements represented by the \ref Intrusive_Link components.
 */
void destructor_IL
	(struct Intrusive_List* lst ///< Standard.
	);

/// \brief Add an \ref Intrusive_Link to the end of the \ref Intrusive_List.
void push_back_IL
	(struct Intrusive_List* lst, ///< The list.
	 struct Intrusive_Link* curr ///< The current link.
	);

/** \brief Erase the current \ref Intrusive_Link.
 *
 *	\warning Does not free memory.
 *
 *	\return Pointer to the next \ref Intrusive_Link.
 */
struct Intrusive_Link* erase_IL
	(struct Intrusive_List* lst, ///< The list.
	 struct Intrusive_Link* curr ///< The current link.
	);

#endif // DPG__Intrusive_h__INCLUDED
