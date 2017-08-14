// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "Matrix.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Macros.h"
#include "allocators.h"

// Will never have _# where # > 1 as this would require storing the dimension of the other levels of dereferencing
// (which should be in the Multiarray). Further, as we consider passing between const/non-const, we want all allocations
// in dynamic memory requiring # > 0. Together these imply that only # = 1 functions are required and the subscript is
// redundant (Same for multiarrays). \todo Add a comment
static struct Matrix_d make_local_Matrix_d_0
	(const char layout, const size_t ext_0, const size_t ext_1, const bool owns_data, double*const data);
static struct Matrix_ui make_local_Matrix_ui_0
	(const char layout, const size_t ext_0, const size_t ext_1, const bool owns_data, unsigned int*const data);

// Constructor/Destructor functions ********************************************************************************* //

struct Matrix_d* constructor_empty_Matrix_d (const char layout, const size_t ext_0, const size_t ext_1)
{
	double* data = mallocator(DOUBLE_T,1,ext_0*ext_1); // keep

	struct Matrix_d local = make_local_Matrix_d_0(layout,ext_0,ext_1,true,data);

	struct Matrix_d* a = malloc(sizeof *a); // returned
	memcpy(a,&local,sizeof *a);

	return a;
}

void const_constructor_move_Matrix_d (const struct const_Matrix_d*const* dest, struct Matrix_d* src)
{
	*(struct const_Matrix_d**) dest = (struct const_Matrix_d*) src;
}

void destructor_Matrix_d (struct Matrix_d* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

struct Matrix_ui* constructor_empty_Matrix_ui (const char layout, const size_t ext_0, const size_t ext_1)
{
	unsigned int* data = mallocator(UINT_T,1,ext_0*ext_1); // keep

	struct Matrix_ui local = make_local_Matrix_ui_0(layout,ext_0,ext_1,true,data);

	struct Matrix_ui* a = malloc(sizeof *a); // returned
	memcpy(a,&local,sizeof *a);

	return a;
}

void const_constructor_move_Matrix_ui (const struct const_Matrix_ui*const* dest, struct Matrix_ui* src)
{
	*(struct const_Matrix_ui**) dest = (struct const_Matrix_ui*) src;
}

void destructor_Matrix_ui (struct Matrix_ui* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

// Helper functions ************************************************************************************************* //

double* get_row_Matrix_d (const size_t row, const struct Matrix_d* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->extents[1])];
}

unsigned int* get_row_Matrix_ui (const size_t row, const struct Matrix_ui* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->extents[1])];
}

// Printing functions *********************************************************************************************** //

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

void print_const_Matrix_d (const struct const_Matrix_d*const a)
{
	struct Matrix_d local = make_local_Matrix_d_0(a->layout,a->extents[0],a->extents[1],false,(double*)a->data);
	print_Matrix_d(&local);
}

void print_Matrix_ui (const struct Matrix_ui*const a)
{
	const size_t ext_r = a->extents[0],
	             ext_c = a->extents[1];

	const unsigned int* data = a->data;

	switch (a->layout) {
	case 'R':
		for (size_t i = 0; i < ext_r; i++) {
			for (size_t j = 0; j < ext_c; j++) {
				const unsigned int val = *data++;
				printf("% 12d ",val);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (size_t i = 0; i < ext_r; i++) {
			for (size_t j = 0; j < ext_c; j++) {
				const unsigned int val = data[i+ext_r*j];
				printf("% 12d ",val);
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

void print_const_Matrix_ui (const struct const_Matrix_ui*const a)
{
	struct Matrix_ui local = make_local_Matrix_ui_0(a->layout,a->extents[0],a->extents[1],false,(unsigned int*)a->data);
	print_Matrix_ui(&local);
}

// Static functions ************************************************************************************************* //

/// \brief Make a local copy of a \ref Matrix_d (static memory).
static struct Matrix_d make_local_Matrix_d_0
	(const char layout,    ///< Standard.
	 const size_t ext_0,   ///< Standard.
	 const size_t ext_1,   ///< Standard.
	 const bool owns_data, ///< Standard.
	 double*const data     ///< Standard.
	)
{
	struct Matrix_d local = { .layout     = layout,
	                          .extents[0] = ext_0,
	                          .extents[1] = ext_1,
	                          .owns_data  = owns_data,
	                          .data       = data,
	                        };
	return local;
}

/// \brief Make a local copy of a \ref Matrix_ui (static memory).
static struct Matrix_ui make_local_Matrix_ui_0
	(const char layout,      ///< Standard.
	 const size_t ext_0,     ///< Standard.
	 const size_t ext_1,     ///< Standard.
	 const bool owns_data,   ///< Standard.
	 unsigned int*const data ///< Standard.
	)
{
	struct Matrix_ui local = { .layout     = layout,
	                           .extents[0] = ext_0,
	                           .extents[1] = ext_1,
	                           .owns_data  = owns_data,
	                           .data       = data,
	                         };
	return local;
}

