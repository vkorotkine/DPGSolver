// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "Intrusive.h"

#include <stddef.h>
#include <stdlib.h>

void init_IL (struct Intrusive_List* lst)
{
	lst->first = NULL;
	lst->last  = NULL;
}

struct Intrusive_List* constructor_empty_IL ()
{
	struct Intrusive_List* lst = malloc(sizeof *lst); // keep
	init_IL(lst);
	return lst;
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
