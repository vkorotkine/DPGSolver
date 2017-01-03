// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__array_swap_h__INCLUDED
#define DPG__array_swap_h__INCLUDED

#include <complex.h>


extern void array_swap_ui         (register unsigned int   *arr1, register unsigned int   *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap_d          (register double         *arr1, register double         *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap_cmplx      (register double complex *arr1, register double complex *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_rearrange_d     (const unsigned int NRows, const unsigned int NCols, const unsigned int *Ordering, const char layout, double *A);
extern void array_rearrange_cmplx (const unsigned int NRows, const unsigned int NCols, const unsigned int *Ordering, const char layout, double complex *A);

#endif // DPG__array_swap_h__INCLUDED
