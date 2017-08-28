// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "allocators.h"

#include <stdlib.h>
#include <stdarg.h>

#include "Macros.h"

static void* mallocator_local (const enum TYPE type, const size_t order, const size_t*const N);
void deallocator_local (void* a, const enum TYPE type, const size_t order, const size_t*const N);

void* mallocator (const enum TYPE type, const size_t order, ...)
{
	if (!order)
		EXIT_UNSUPPORTED;

	va_list ap;
	va_start(ap,order); // free

	size_t N[order];
	for (size_t i = 0; i < order; i++)
		N[i] = va_arg(ap,size_t);

	va_end(ap);

	return mallocator_local(type,order,N);
}

void deallocator (void* a, const enum TYPE type, const size_t order, ...)
{
	if (!order)
		EXIT_UNSUPPORTED;

	va_list ap;
	va_start(ap,order); // free

	size_t N[order];
	for (size_t i = 0; i < order; i++)
		N[i] = va_arg(ap,size_t);

	va_end(ap);

	deallocator_local(a,type,order,N);
}


// Static functions

static void* mallocator_local (const enum TYPE type, const size_t order, const size_t*const N)
{
	// It was not possible to make this function more compact as the 'variable type' cannot be passed.
	const size_t index = order-1;

	switch (type) {
	case CHAR_T:
		switch (order) {
		case 1: {
			char* A = malloc(N[index] * sizeof *A);
			return A;
			break;
		} case 2: {
			char** A = malloc(N[index] * sizeof *A);
			for (size_t i = 0; i < N[index]; i++)
				A[i] = mallocator_local(CHAR_T,index,N);
			return A;
			break;
		} case 0: default:
			EXIT_UNSUPPORTED;
			break;
		}
		break;
	case UINT_T:
		switch (order) {
		case 1: {
			unsigned int* A = malloc(N[index] * sizeof *A);
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
	case SIZE_T_T:
		switch (order) {
		case 1: {
			size_t* A = malloc(N[index] * sizeof *A);
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

void deallocator_local (void* a, const enum TYPE type, const size_t order, const size_t*const N)
{
	const size_t index = order-1;

	switch (type) {
	case CHAR_T:
		switch (order) {
		case 1: {
			char* b = a;
			free(b);
			break;
		} case 2: {
			char** b = a;
			for (size_t i = 0; i < N[index]; i++)
				deallocator_local(b[i],type,index,N);
			free(b);
			break;
		} case 0:
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
