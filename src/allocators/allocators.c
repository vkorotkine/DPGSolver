// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "allocators.h"

#include <stdlib.h>
#include <stdarg.h>

#include "Macros.h"

// Static function declarations ************************************************************************************* //

/// \brief Perform allocation of arrays of standard data types using malloc.
static void* mallocator_local
	(const enum Variable_Type type, ///< Defined in \ref mallocator.
	 const int order,               ///< Defined in \ref mallocator.
	 const ptrdiff_t*const N        ///< The array of variadic arguments passed to \ref mallocator.
	);

/// \brief Perform deallocation of arrays of standard data types.
static void deallocator_local
	(void* a,                       ///< Define in \ref deallocator.
	 const enum Variable_Type type, ///< Define in \ref deallocator.
	 const int order,               ///< Define in \ref deallocator.
	 const ptrdiff_t*const N        ///< The array of variadic arguments passed to \ref deallocator.
	);

// Interface functions ********************************************************************************************** //

void* mallocator (const enum Variable_Type type, const int order, ...)
{
	if (!order)
		EXIT_UNSUPPORTED;

	va_list ap;
	va_start(ap,order); // free

	ptrdiff_t N[order];
	for (ptrdiff_t i = 0; i < order; i++)
		N[i] = va_arg(ap,ptrdiff_t);

	va_end(ap);

	return mallocator_local(type,order,N);
}

void deallocator (void* a, const enum Variable_Type type, const int order, ...)
{
	if (!order)
		EXIT_UNSUPPORTED;

	va_list ap;
	va_start(ap,order); // free

	ptrdiff_t N[order];
	for (ptrdiff_t i = 0; i < order; i++)
		N[i] = va_arg(ap,ptrdiff_t);

	va_end(ap);

	deallocator_local(a,type,order,N);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void* mallocator_local (const enum Variable_Type type, const int order, const ptrdiff_t*const N)
{
	const ptrdiff_t index = order-1;

	switch (type) {
	case CHAR_T:
		switch (order) {
		case 1: {
			char* A = malloc(N[index] * sizeof *A);
			return A;
			break;
		} case 2: {
			char** A = malloc(N[index] * sizeof *A);
			for (ptrdiff_t i = 0; i < N[index]; i++)
				A[i] = mallocator_local(type,index,N);
			return A;
			break;
		} case 0: default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case INT_T:
		switch (order) {
		case 1: {
			int* A = malloc(N[index] * sizeof *A);
			return A;
			break;
		} case 0: default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case DOUBLE_T:
		switch (order) {
		case 1: {
			double* A = malloc(N[index] * sizeof *A);
			return A;
			break;
		} case 0: default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case PTRDIFF_T:
		switch (order) {
		case 1: {
			ptrdiff_t* A = malloc(N[index] * sizeof *A);
			return A;
			break;
		} case 0: default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

static void deallocator_local (void* a, const enum Variable_Type type, const int order, const ptrdiff_t*const N)
{
	const ptrdiff_t index = order-1;

	if (order == 1) {
		free(a);
		return;
	}

	switch (type) {
	case CHAR_T:
		switch (order) {
		case 2: {
			char** b = a;
			for (ptrdiff_t i = 0; i < N[index]; i++)
				deallocator_local(b[i],type,index,N);
			free(b);
			break;
		} case 0: case 1:
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case INT_T:
		switch (order) {
		case 0: case 1:
		default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}
