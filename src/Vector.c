// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "Vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

struct Vector_ui* constructor_empty_Vector_ui (const size_t n_rows)
{
	unsigned int* data = malloc(n_rows * sizeof *data); // keep

	struct Vector_ui local = { .extents[0] = n_rows,
	                           .owns_data  = true,
	                           .data       = data,
	                         };

	struct Vector_ui* a = malloc(sizeof *a); // returned
	memcpy(a,&local,sizeof *a);

	return a;
}

void const_constructor_const_Vector_ui_2_Vector_ui
	(const struct Vector_ui*const*const* dest, struct Vector_ui** src, const size_t n0)
{
	// See Matrix.c for implementation outline.
}

void reorder_Vector_ui (struct Vector_ui*const a, unsigned int*const ordering)
{
	const size_t size = a->extents[0];

	unsigned int b[size];
	for (size_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (size_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void print_Vector_ui (const struct Vector_ui*const a)
{
	const size_t ext = a->extents[0];

	const unsigned int* data = a->data;

	for (size_t i = 0; i < ext; i++)
		printf("% 12d ",*data++);
	printf("\n");
}
