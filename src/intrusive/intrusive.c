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
/// \file

#include "intrusive.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_intrusive.h"

// Static function declarations ************************************************************************************* //

/** \brief Return the pointer to the link for which \ref Intrusive_Link::next == NULL.
 *  \return See brief. */
static struct Intrusive_Link* get_last_link
	(struct Intrusive_Link*const first ///< The first link.
	);

/** \brief Return the pointer to the link `n` entries after the input.
 *  \return See brief. */
static struct Intrusive_Link* advance_link
	(struct Intrusive_Link* curr, ///< Pointer to the current link.
	 const int n                  ///< The number of links to skip over.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_empty_IL (const int list_name, struct Intrusive_List* base)
{
	struct Intrusive_List* lst = calloc(1,sizeof *lst); // keep

	lst->first = NULL;
	lst->last  = NULL;

	lst->name  = list_name;
	lst->base  = base;

	return lst;
}

const struct const_Intrusive_List* constructor_empty_const_IL
	(const int list_name, const struct const_Intrusive_List* base)
{
	return (const struct const_Intrusive_List*) constructor_empty_IL(list_name,(struct Intrusive_List*)base);
}

void clear_IL (struct Intrusive_List* lst)
{
	struct Intrusive_Link* curr = lst->first;
	while (curr) {
		struct Intrusive_Link* next = curr->next;
		free(curr);
		curr = next;
	}
	lst->first = NULL;
	lst->last  = NULL;
}

void destructor_IL (struct Intrusive_List* lst, const bool destruct_links)
{
	if (destruct_links)
		clear_IL(lst);
	free(lst);
}

void destructor_const_IL (const struct const_Intrusive_List* lst, const bool destruct_links)
{
	destructor_IL((struct Intrusive_List*)lst,destruct_links);
}

void destructor_IL_base (struct Intrusive_List* lst)
{
	assert(lst->base);
	destructor_IL(lst->base,true);

	lst->base = NULL;
}

void destructor_const_IL_base (const struct const_Intrusive_List* lst)
{
	destructor_IL_base((struct Intrusive_List*)lst);
}

void push_back_IL (struct Intrusive_List* lst, struct Intrusive_Link* curr)
{
	struct Intrusive_Link* last = lst->last;
	if (last) {
		// Add `curr` after last.
		last->next = curr;
		curr->prev = last;
	} else {
		// `curr` is the first element.
		lst->first = curr;
		curr->prev = NULL;
	}
	// `curr` is the last element.
	lst->last = curr;
	curr->next = NULL;
}

void push_back_const_IL (const struct const_Intrusive_List* lst, const struct const_Intrusive_Link* curr)
{
	push_back_IL((struct Intrusive_List*)lst,(struct Intrusive_Link*)curr);
}

struct Intrusive_Link* erase_IL (struct Intrusive_List* lst, struct Intrusive_Link* curr)
{
	if (curr == lst->first) {
		if (curr->next) {
			lst->first = curr->next;
			curr->next->prev = NULL;
			return curr->next;
		} else {
			lst->first = NULL;
			lst->last  = NULL;
			return NULL;
		}
	} else if (curr == lst->last) {
		if (curr->prev) {
			lst->last = curr->prev;
			curr->prev->next = NULL;
			return NULL;
		} else {
			lst->first = NULL;
			lst->last  = NULL;
			return NULL;
		}
	} else {
		curr->next->prev = curr->prev;
		curr->prev->next = curr->next;
		return curr->next;
	}
}

const struct const_Intrusive_Link* erase_const_IL
	(const struct const_Intrusive_List* lst, const struct const_Intrusive_Link* curr)
{
	return (const struct const_Intrusive_Link*) erase_IL((struct Intrusive_List*)lst,(struct Intrusive_Link*)curr);
}

struct Intrusive_Link* replace_IL
	(struct Intrusive_List*const lst, const int n_replace, struct Intrusive_Link*const curr_f,
	 struct Intrusive_Link*const first)
{
	struct Intrusive_Link*const last   = get_last_link(first);
	struct Intrusive_Link*const curr_l = advance_link(curr_f,n_replace-1);
	if (curr_f == lst->first) {
		lst->first = first;
		first->prev = NULL;
		last->next = curr_l->next;
		if (curr_l->next) {
			curr_l->next->prev = last;
		} else {
			lst->last = last;
		}
	} else if (curr_l == lst->last) {
		lst->last  = last;
		assert(last->next == NULL);
		first->prev = curr_f->prev;
		if (curr_f->prev) {
			curr_f->prev->next = first;
		} else {
			lst->first = first;
		}
	} else {
		last->next  = curr_l->next;
		curr_l->next->prev = last;

		first->prev = curr_f->prev;
		curr_f->prev->next = first;
	}
	return last;
}

struct Intrusive_Link* advance_Link (const int n_adv, struct Intrusive_Link* first)
{
	struct Intrusive_Link* curr = first;
	for (int n = 0; n < n_adv; ++n)
		curr = curr->next;
	return curr;
}

void insert_List_into_List
	(struct Intrusive_List*const sub, struct Intrusive_List*const main, struct Intrusive_Link*const curr)
{
	if (curr == main->first) {
		sub->first->prev = NULL;
		main->first = sub->first;
	} else {
		sub->first->prev = curr->prev;
		curr->prev->next = sub->first;
	}
	sub->last->next = curr;
	curr->prev = sub->last;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Intrusive_Link* get_last_link (struct Intrusive_Link*const first)
{
	struct Intrusive_Link* curr = first;
	while (curr->next)
		curr = curr->next;
	return curr;
}

static struct Intrusive_Link* advance_link (struct Intrusive_Link* curr, const int n)
{
	for (int i = 0; i < n; ++i) {
		assert(curr->next != NULL);
		curr = curr->next;
	}
	return curr;
}
