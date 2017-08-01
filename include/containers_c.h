// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__containers_c_h__INCLUDED
#define DPG__containers_c_h__INCLUDED
/**	\file
 *	\brief Provides containers and related functions for `complex` datatypes.
 *
 *	\see Multiarray_d provides the documentation.
 */

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

///	\brief Supports dense multi-dimensional (double complex) arrays.
struct Multiarray_c {
	char layout;

	size_t  order;
	size_t* extents;

	bool            owns_data;
	double complex* data;
};

struct Multiarray_c* constructor_empty_Multiarray_c_1 (const char layout, const size_t order, ...);

void destructor_Multiarray_c_1 (struct Multiarray_c* A);

#endif // DPG__containers_c_h__INCLUDED
