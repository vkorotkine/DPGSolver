// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "containers.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "Macros.h"

// Multiarray_* ***************************************************************************************************** //

struct Multiarray_d* constructor_move_Multiarray_d_1_d (const char layout, double*const data, const size_t order, ...)
{
	va_list ap;
	va_start(ap,order); // free
	size_t* extents = set_extents(order,ap); // keep
	va_end(ap);

	struct Multiarray_d local = { .layout    = layout,
	                              .order     = order,
	                              .extents   = extents,
	                              .owns_data = false,
	                              .data      = data,
	                            };

	struct Multiarray_d* A = malloc(sizeof *A); // returned
	memcpy(A,&local,sizeof *A);

	return A;
}

void destructor_Multiarray_d_1 (struct Multiarray_d* A)
{
	free(A->extents);
	if (A->owns_data)
		free(A->data);
	FREE_NULL(A);
}

size_t* set_extents (const size_t order, va_list ap)
{
	size_t* extents = malloc(order * sizeof *extents); // returned

	for (size_t i = 0; i < order; i++)
		extents[i] = va_arg(ap,size_t);

	return extents;
}

size_t compute_size (const size_t order, const size_t *const extents)
{
	size_t size = 1;
	for (size_t i = 0; i < order; i++)
		size *= extents[i];

	return size;
}

// Matrix_* ********************************************************************************************************* //

struct Matrix_d* constructor_empty_Matrix_d (const char layout, const size_t n_rows, const size_t n_cols)
{
	double* data = malloc(n_rows*n_cols * sizeof *data); // keep

	struct Matrix_d local = { .layout     = layout,
	                          .extents[0] = n_rows,
	                          .extents[1] = n_cols,
	                          .owns_data  = true,
	                          .data       = data,
	                        };

	struct Matrix_d* a = malloc(sizeof *a); // returned
	memcpy(a,&local,sizeof *a);

	return a;
}

double* get_row_Matrix_d (const size_t row, const struct Matrix_d* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->extents[1])];
}

void print_Matrix_d (const struct Matrix_d*const a)
{
	const size_t ext_r = a->extents[0],
	             ext_c = a->extents[1];

	const double* data = a->data;

	switch (a->layout) {
	case 'R':
		for (size_t i = 0; i < ext_r; i++) {
			for (size_t j = 0; j < ext_c; j++) {
				const double val = *data++;
				printf("% .4e ",( (isnan(val) || fabs(val)) ? val : 0.0 ));
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (size_t i = 0; i < ext_r; i++) {
			for (size_t j = 0; j < ext_c; j++) {
				const double val = data[i+ext_r*j];
				printf("% .4e ",( (isnan(val) || fabs(val)) ? val : 0.0 ));
			}
			printf("\n");
		}
		printf("\n");
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

// Vector_* ********************************************************************************************************* //

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

void print_Vector_ui (const struct Vector_ui*const a)
{
	const size_t ext = a->extents[0];

	const unsigned int* data = a->data;

	for (size_t i = 0; i < ext; i++)
		printf("% 12d ",*data++);
	printf("\n");
}
