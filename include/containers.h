// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__containers_h__INCLUDED
#define DPG__containers_h__INCLUDED

#include <stddef.h>
#include <stdarg.h>

/*	Constructors for const_Multiarray_* structs are not provided directly. To assign to these variables, an lvalue cast
 *	to the appropriate non-const type should be used. This procedure does not exhibit undefined behaviour (relating to
 *	changing a const-qualified type) as long as the memory for the structs is dynamically allocated.
 *
 *	Example:
 *		To assign to `const struct const_Multiarray_d*const var1`, use the relevant constructor to set a variable:
 *		`struct Multiarray_d* var2`, then cast the lvalue `*(struct Multiarray_d**)&var1 = var2`.
 */

struct Multiarray_d {
	/*	Purpose:
	 *		Defines a struct to support dense multi-dimensional (double) arrays.
	 *
	 *	Comments:
	 *		The Multi_array struct is intended to be used as a higher-dimensional matrix where move constructors are
	 *		used to form matrix structs for appropriate sub-blocks. As the data is stored contiguously in memory, the
	 *		multi-array may also be acted on over multiple dimensions at once.
	 */

	char layout; // May be (R)ow or (C)olumn major.

	size_t  order;   // Number of dimensions.
	size_t* extents; // Size of arrays in each dimension.

	double* data; // The data.
};

struct const_Multiarray_d {
	/*	Purpose:
	 *		const version of Multiarray_d.
	 */

	const char layout;

	const size_t       order;
	const size_t*const extents;

	const double*const data;
};

struct Multiarray_d* constructor_move_Multiarray_d_1_d (const char layout, double*const data, const size_t order, ...);

void destructor_moved_Multiarray_d_1 (struct Multiarray_d* A);


size_t* set_extents  (const size_t order, va_list ap);
size_t  compute_size (const size_t n_ext, const size_t *const extents);


#endif // DPG__containers_h__INCLUDED
