// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "vector.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"

#include "allocators.h"
#include "multiarray.h"
#include "matrix.h"

// Static function declarations ************************************************************************************* //

/** \brief Make a local \ref Vector_i\* (dynamic memory).
 *	\return See brief.  */
static struct Vector_i* constructor_local_Vector_i_1
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 int*const data         ///< Standard.
	);

/** \brief Make a local \ref Vector_d\* (dynamic memory).
 *	\return See brief.  */
static struct Vector_d* constructor_local_Vector_d_1
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 double*const data      ///< Standard.
	);

/** \brief Comparison function for std::qsort between `int*` `a` and `b`.
 *	\return a - b.  */
static int cmp_i
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/** \brief Constructs a default \ref Vector_i\*.
 *	\return Standard. */
static struct Vector_i* constructor_default_Vector_i ();

// Constructor/Destructor functions ********************************************************************************* //

struct Vector_i** constructor_default_Vector_i_2 (const ptrdiff_t n_dest)
{
	struct Vector_i** dest = malloc(n_dest * sizeof *dest); // returned;

	for (ptrdiff_t n = 0; n < n_dest; n++)
		dest[n] = constructor_default_Vector_i();

	return dest;
}

struct Vector_i* constructor_empty_Vector_i (const ptrdiff_t ext_0)
{
	int* data = mallocator(INT_T,1,ext_0); // keep

	return constructor_local_Vector_i_1(ext_0,true,data);
}

struct Vector_i* constructor_copy_Vector_i (const struct Vector_i*const src)
{
	const ptrdiff_t ext_0 = src->ext_0;
	const int*const data_src = src->data;

	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_local_Vector_i_1(ext_0,true,data);
}

struct Vector_i* constructor_copy_Vector_i_i (const ptrdiff_t ext_0, const int*const data_src)
{
	int* data = malloc(ext_0 * sizeof *data); // keep
	for (ptrdiff_t i = 0; i < ext_0; i++)
		data[i] = data_src[i];

	return constructor_local_Vector_i_1(ext_0,true,data);
}

struct Vector_i* constructor_move_Vector_i_i (const ptrdiff_t ext_0, const bool owns_data, int*const data)
{
	return constructor_local_Vector_i_1(ext_0,owns_data,data);
}

struct const_Vector_i* constructor_move_const_Vector_i_i
	(const ptrdiff_t ext_0, const bool owns_data, const int*const data)
{
	struct Vector_i* local = constructor_local_Vector_i_1(ext_0,owns_data,(int*)data); // free

	struct const_Vector_i* dest = malloc(sizeof *dest); // returned
	memcpy(dest,local,sizeof *dest);
	free(local);

	return dest;
}

void const_constructor_move_Vector_i (const struct const_Vector_i*const* dest, struct Vector_i* src)
{
	*(struct const_Vector_i**) dest = (struct const_Vector_i*) src;
}

void destructor_Vector_i (struct Vector_i* a)
{
	if (a->owns_data)
		deallocator(a->data,INT_T,1,a->ext_0);
	free(a);
}

void destructor_Vector_i_2 (struct Vector_i** a, const ptrdiff_t n_src, const bool owns_data)
{
	if (owns_data) {
		for (ptrdiff_t n = 0; n < n_src; n++)
			destructor_Vector_i(a[n]);
	}
	free(a);
}

struct Vector_d* constructor_empty_Vector_d (const ptrdiff_t ext_0)
{
	double* data = mallocator(DOUBLE_T,1,ext_0); // keep

	return constructor_local_Vector_d_1(ext_0,true,data);
}

struct Vector_d* constructor_sum_Vector_d_const_Matrix_d (const char sum_dir, const struct const_Matrix_d*const src)
{
	if (!(sum_dir == 'R' || sum_dir == 'C'))
		EXIT_UNSUPPORTED;

	const ptrdiff_t ext_0 = ( sum_dir == 'R' ? src->ext_1 : src->ext_0 );

	struct Vector_d* dest = constructor_empty_Vector_d(ext_0); // returned
	for (ptrdiff_t j = 0; j < ext_0; ++j)
		dest->data[j] = 0.0;

	if (sum_dir != src->layout) {
		EXIT_ADD_SUPPORT;
	} else {
		if (src->layout == 'R') {
			const ptrdiff_t i_max = src->ext_0;
			for (ptrdiff_t i = 0; i < i_max; ++i) {
				const double* data_m = get_row_const_Matrix_d(i,src);
				for (ptrdiff_t j = 0; j < ext_0; ++j)
					dest->data[j] += data_m[j];
			}
		} else {
			EXIT_ADD_SUPPORT;
		}
	}

	return dest;
}

void destructor_Vector_d (struct Vector_d* a)
{
	if (a->owns_data)
		deallocator(a->data,DOUBLE_T,1,a->ext_0);
	free(a);
}

// Helper functions ************************************************************************************************* //

void reorder_Vector_i (struct Vector_i*const a, const int*const ordering)
{
	const ptrdiff_t size = a->ext_0;

	int b[size];
	for (ptrdiff_t i = 0; i < size; i++)
		b[i] = a->data[ordering[i]];

	for (ptrdiff_t i = 0; i < size; i++)
		a->data[i] = b[i];
}

void resize_Vector_i (struct Vector_i*const a, const ptrdiff_t ext_0)
{
	const ptrdiff_t size_i = a->ext_0;
	a->ext_0 = ext_0;
	const ptrdiff_t size_o = a->ext_0;

	if (size_o <= size_i)
		return;

	const int* data_i = a->data;
	a->data = malloc(size_o * sizeof *(a->data)); // keep
	if (size_i != 0) {
		for (ptrdiff_t i = 0; i < size_i; i++)
			a->data[i] = data_i[i];
	}
	free((void*)data_i);
}

void set_to_zero_Vector_i (struct Vector_i*const a)
{
	const ptrdiff_t i_max = a->ext_0;
	for (ptrdiff_t i = 0; i < i_max; i++)
		a->data[i] = 0;
}

void set_to_data_Vector_i (struct Vector_i*const a, const int*const data_src)
{
	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		a->data[i] = data_src[i];
}

void sort_Vector_i (struct Vector_i* a)
{
	const ptrdiff_t size = a->ext_0;
	qsort(a->data,size,sizeof(a->data[0]),cmp_i);
}

int sum_Vector_i (struct Vector_i* a)
{
	int sum = 0;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; ++i)
		sum += a->data[i];
	return sum;
}

bool check_equal_Vector_i (const struct Vector_i*const a, const struct Vector_i*const b)
{
	const ptrdiff_t size = a->ext_0;
	if (size != b->ext_0)
		return false;

	const int* data_a = a->data,
	         * data_b = b->data;

	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

bool check_equal_Vector_i_i (const struct Vector_i*const a, const int* data_b)
{
	const int* data_a = a->data;

	const ptrdiff_t size = a->ext_0;
	for (ptrdiff_t i = 0; i < size; i++) {
		if (*data_a++ != *data_b++)
			return false;
	}
	return true;
}

int cmp_Vector_i (const void *a, const void *b)
{
	const struct Vector_i*const*const ia = (const struct Vector_i*const*const) a,
	                     *const*const ib = (const struct Vector_i*const*const) b;

	const ptrdiff_t size_a = (*ia)->ext_0,
	                size_b = (*ib)->ext_0;

	if (size_a > size_b)
		return 1;
	else if (size_a < size_b)
		return -1;

	const int*const data_a = (*ia)->data,
	         *const data_b = (*ib)->data;

	for (ptrdiff_t i = 0; i < size_a; ++i) {
		if (data_a[i] > data_b[i])
			return 1;
		else if (data_a[i] < data_b[i])
			return -1;
	}
	return 0;
}

void copy_data_Vector_i_Vector_i (const struct Vector_i*const src, struct Vector_i*const dest)
{
	const ptrdiff_t size_src  = src->ext_0,
	                size_dest = dest->ext_0;

	if (size_src != size_dest)
		EXIT_UNSUPPORTED;

	for (ptrdiff_t i = 0; i < size_src; ++i)
		dest->data[i] = src->data[i];
}

void push_back_Vector_i (struct Vector_i*const src, const int val, const bool sorted, const bool unique)
{
	if (sorted)
		sort_Vector_i(src);

	const bool add_val = unique ? !find_val_Vector_i((struct const_Vector_i*)src,val,sorted) : true;

	if (!add_val)
		return;

	resize_Vector_i(src,src->ext_0+1);
	src->data[src->ext_0-1] = val;

	if (sorted)
		sort_Vector_i(src);
}

bool find_val_Vector_i (const struct const_Vector_i*const src, const int val, const bool sorted)
{
	bool found = false;
	if (!sorted) {
		const ptrdiff_t i_max = src->ext_0;
		for (ptrdiff_t i = 0; i < i_max; ++i) {
			if (src->data[i] == val) {
				found = true;
				break;
			}
		}
	} else {
		const int* ind_ptr = bsearch(&val,src->data,src->ext_0,sizeof(src->data[0]),cmp_i);
		if (ind_ptr)
			found = true;
	}
	return found;
}

// Printing functions *********************************************************************************************** //

void print_Vector_i (const struct Vector_i*const a)
{
	const ptrdiff_t ext = a->ext_0;

	const int* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
		printf("% 12d ",*data++);
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

void print_const_Vector_i (const struct const_Vector_i*const a)
{
	struct Vector_i* local = constructor_local_Vector_i_1(a->ext_0,false,(int*)a->data); // free
	print_Vector_i(local);
	free(local);
}

void print_Vector_d (const struct Vector_d*const a, const double tol)
{
	const ptrdiff_t ext = a->ext_0;

	const double* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
		const double val = *data++;
		printf("% .4e ",( (isnan(val) || (fabs(val) > tol)) ? val : 0.0 ));
		if (!((i+1)%8))
			printf("\n");
	}
	printf("\n\n");
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static struct Vector_i* constructor_default_Vector_i ()
{
	struct Vector_i* dest = malloc(sizeof *dest); // returned
	dest->ext_0     = 0;
	dest->owns_data = true;
	dest->data      = NULL;

	return dest;
}

static struct Vector_i* constructor_local_Vector_i_1 (const ptrdiff_t ext_0, const bool owns_data, int*const data)
{
	struct Vector_i* dest = malloc(sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

static struct Vector_d* constructor_local_Vector_d_1 (const ptrdiff_t ext_0, const bool owns_data, double*const data)
{
	struct Vector_d* dest = malloc(sizeof *dest); // returned

	dest->ext_0     = ext_0;
	dest->owns_data = owns_data;
	dest->data      = data;

	return dest;
}

static int cmp_i (const void *a, const void *b)
{
	return (int) ( *(int*)a - *(int*)b );
}
