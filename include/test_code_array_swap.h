// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_array_swap_h__INCLUDED
#define DPG__test_code_array_swap_h__INCLUDED

extern void array_swap1_ui (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap2_ui (register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap1_d  (register double       *arr1, register double       *arr2, const unsigned int NIn, const unsigned int stepIn);
extern void array_swap2_d  (register double       *arr1, register double       *arr2, const unsigned int NIn, const unsigned int stepIn);

#endif // DPG__test_code_array_swap_h__INCLUDED
