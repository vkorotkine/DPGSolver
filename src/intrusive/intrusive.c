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

void destructor_IL (struct Intrusive_List* lst)
{
	clear_IL(lst);
	free(lst);
}

void destructor_const_IL (const struct const_Intrusive_List* lst)
{
	destructor_IL((struct Intrusive_List*)lst);
}

void destructor_IL_base (struct Intrusive_List* lst)
{
	assert(lst->base);
	destructor_IL(lst->base);

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

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
