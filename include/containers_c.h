// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__containers_c_h__INCLUDED
#define DPG__containers_c_h__INCLUDED

#include <stddef.h>
#include <complex.h>

struct Multiarray_c {
	/*
	 *	Purpose:
	 *		Defines a struct to support dense multi-dimensional (double complex) arrays.
	 *
	 *	Comments:
	 *		The Multi_array struct is intended to be used as a higher-dimensional matrix where move constructors are
	 *		used to form matrix structs for appropriate sub-blocks. As the data is stored contiguously in memory, the
	 *		multi-array may also be acted on over multiple dimensions at once.
	 */

	char layout; // May be (R)ow or (C)olumn major.

	size_t  order;   // Number of dimensions.
	size_t* extents; // Size of arrays in each dimension

	double complex* data;
};

struct Multiarray_c* constructor_empty_Multiarray_c_1 (const char layout, const size_t order, ...);

#endif // DPG__containers_c_h__INCLUDED
