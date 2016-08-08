// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__array_swap_h__INCLUDED
#define DPG__array_swap_h__INCLUDED

#include <stdlib.h>
 

extern void array_swap_ui (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap_d  (register double *arr1,       register double *arr2,       const int NIn,          const int stepIn);

#endif // DPG__array_swap_h__INCLUDED
