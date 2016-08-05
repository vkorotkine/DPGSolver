// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__array_print_h__INCLUDED
#define DPG__array_print_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>

extern void array_print_ui (const unsigned int m, const unsigned int n, const unsigned int *A, const char layout);
extern void array_print_i  (const unsigned int m, const unsigned int n, int *A,                const char layout);
extern void array_print_l  (const unsigned int m, const unsigned int n, long *A,               const char layout);
extern void array_print_ll (const unsigned int m, const unsigned int n, long long *A,          const char layout);
extern void array_print_f  (const unsigned int m, const unsigned int n, float *A,              const char layout);
extern void array_print_d  (const unsigned int m, const unsigned int n, double *A,             const char layout);
extern void array_print_ld (const unsigned int m, const unsigned int n, long double *A,        const char layout);

#endif // DPG__array_print_h__INCLUDED
