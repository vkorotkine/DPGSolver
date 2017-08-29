// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/** \file
 */

#include "matrix.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Macros.h"
#include "allocators.h"

// Static function declarations ************************************************************************************* //

/** \brief Make a local \ref Matrix_d\* (dynamic memory).
 *	\return See brief. */
static struct Matrix_d* constructor_local_Matrix_d_1
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/** \brief Make a local \ref Matrix_i\* (dynamic memory).
 *	\return See brief. */
static struct Matrix_i* constructor_local_Matrix_i_1
	(const char layout,     ///< Standard.
	 const ptrdiff_t ext_0, ///< Standard.
	 const ptrdiff_t ext_1, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

// Constructor/Destructor functions ********************************************************************************* //

struct Matrix_d* constructor_empty_Matrix_d (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	double* data = mallocator(DOUBLE_T,1,ext_0*ext_1); // keep
	struct Matrix_d* a = constructor_local_Matrix_d_1(layout,ext_0,ext_1,true,data); // returned

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

struct Matrix_i* constructor_empty_Matrix_i (const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	int* data = mallocator(INT_T,1,ext_0*ext_1); // keep
	struct Matrix_i* a = constructor_local_Matrix_i_1(layout,ext_0,ext_1,true,data); // returned

	return a;
}

void const_constructor_move_Matrix_i (const struct const_Matrix_i*const* dest, struct Matrix_i* src)
{
	*(struct const_Matrix_i**) dest = (struct const_Matrix_i*) src;
}

void destructor_Matrix_i (struct Matrix_i* a)
{
	if (a->owns_data)
		free(a->data);
	free(a);
}

// Helper functions ************************************************************************************************* //

double* get_row_Matrix_d (const ptrdiff_t row, const struct Matrix_d* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->ext_1)];
}

const double* get_row_const_Matrix_d (const ptrdiff_t row, const struct const_Matrix_d*const a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->ext_1)];
}

int* get_row_Matrix_i (const ptrdiff_t row, const struct Matrix_i* a)
{
	if (a->layout != 'R')
		EXIT_UNSUPPORTED;

	return &a->data[row*(a->ext_1)];
}

double compute_norm_Matrix_d_row
	(const ptrdiff_t row, const struct Matrix_d*const a, const char*const norm_type)
{
	const double*const data = get_row_Matrix_d(row,a);
	const ptrdiff_t i_max   = a->ext_1;

	double norm = 0.0;
	if (strstr(norm_type,"L2")) {
		for (ptrdiff_t i = 0; i < i_max; ++i)
			norm += data[i]*data[i];
		return sqrt(norm);
	}
	EXIT_UNSUPPORTED;
}

int get_val_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct Matrix_i*const a)
{
	int*const data = a->data;
	switch (a->layout) {
	case 'R': {
		const ptrdiff_t ext_1 = a->ext_1;
		return data[row*ext_1+col];
	} case 'C':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

int get_val_const_Matrix_i (const ptrdiff_t row, const ptrdiff_t col, const struct const_Matrix_i*const a)
{
	const int*const data = a->data;
	switch (a->layout) {
	case 'R': {
		const ptrdiff_t ext_1 = a->ext_1;
		return data[row*ext_1+col];
	} case 'C':
		EXIT_ADD_SUPPORT;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void set_row_Matrix_d (const ptrdiff_t row, const struct Matrix_d* dest, const double*const data_src)
{
	if (dest->layout != 'R')
		EXIT_ADD_SUPPORT;

	double*const data = get_row_Matrix_d(row,dest);

	const ptrdiff_t i_max = dest->ext_1;
	for (ptrdiff_t i = 0; i < i_max; ++i)
		data[i] = data_src[i];
}

// Printing functions *********************************************************************************************** //

void print_Matrix_d (const struct Matrix_d*const a)
{
	const ptrdiff_t ext_r = a->ext_0,
	                ext_c = a->ext_1;

	const double* data = a->data;

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_r; i++) {
			for (ptrdiff_t j = 0; j < ext_c; j++) {
				const double val = *data++;
				printf("% .4e ",( (isnan(val) || fabs(val)) ? val : 0.0 ));
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_r; i++) {
			for (ptrdiff_t j = 0; j < ext_c; j++) {
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
	struct Matrix_d* local = constructor_local_Matrix_d_1(a->layout,a->ext_0,a->ext_1,false,(double*)a->data);
	print_Matrix_d(local);
	free(local);
}

void print_Matrix_i (const struct Matrix_i*const a)
{
	const ptrdiff_t ext_r = a->ext_0,
	                ext_c = a->ext_1;

	const int* data = a->data;

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_r; i++) {
			for (ptrdiff_t j = 0; j < ext_c; j++) {
				const int val = *data++;
				printf("% 12d ",val);
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_r; i++) {
			for (ptrdiff_t j = 0; j < ext_c; j++) {
				const int val = data[i+ext_r*j];
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

void print_const_Matrix_i (const struct const_Matrix_i*const a)
{
	struct Matrix_i* local = constructor_local_Matrix_i_1(a->layout,a->ext_0,a->ext_1,false,(int*)a->data);
	print_Matrix_i(local);
	free(local);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Matrix_d* constructor_local_Matrix_d_1
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, double*const data)
{
	struct Matrix_d* dest = malloc(sizeof *dest); // returned

	dest->layout    = layout;
	dest->ext_0     = ext_0;
	dest->ext_1     = ext_1;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

static struct Matrix_i* constructor_local_Matrix_i_1
	(const char layout, const ptrdiff_t ext_0, const ptrdiff_t ext_1, const bool owns_data, int*const data)
{
	struct Matrix_i* dest = malloc(sizeof *dest); // returned

	dest->layout    = layout;
	dest->ext_0     = ext_0;
	dest->ext_1     = ext_1;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}
